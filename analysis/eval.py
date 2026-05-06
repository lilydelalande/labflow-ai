"""
eval.py — compare a measurement run against the reference benchmarks.

Two checks (per the design banked in LAB_NOTEBOOK.md):

  Per-image quality check.
      For each new image, compare its dimensionless quality metrics
      (wall_fit_success_rate, reliable_rate, drop_rate, median_wall_cv)
      and capsid_mean_nm against the distribution of the same metrics in
      `benchmarks/<sample_type>/reference_runs.csv`, restricted to the same
      sample_subtype. Flag values outside (median ± 2 × IQR).

  Hand vs script.
      If `benchmarks/<sample_type>/reference_hand.csv` has any rows tagged
      with the new run's batch_id (or — fallback — the same sample_subtype
      as the new run), compute hand_mean − script_mean. This surfaces
      systematic bias between the script and human ImageJ measurements.

Outputs:
  - structured dict (returned)
  - `<run_dir>/eval_report.md` (human-readable)

Usage:
  from analysis.eval import evaluate
  ev = evaluate(run_result_dict_or_run_dir, sample_type="VLP")

  # or CLI:
  uv run python -m analysis.eval results/some_batch
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

REPO_ROOT       = Path(__file__).resolve().parent.parent
DEFAULT_BENCH   = REPO_ROOT / "benchmarks"
IQR_TOL         = 2.0    # flag if metric is more than 2× IQR from median
MIN_N_FOR_FLAG  = 4      # need at least this many ref points to flag (otherwise just report)


# ── Helpers ──────────────────────────────────────────────────────────────

def _load_reference(sample_type: str, benchmarks_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    runs_path = benchmarks_dir / sample_type.lower() / "reference_runs.csv"
    hand_path = benchmarks_dir / sample_type.lower() / "reference_hand.csv"
    if not runs_path.exists():
        raise FileNotFoundError(f"reference_runs.csv not found: {runs_path}")
    runs = pd.read_csv(runs_path)
    hand = pd.read_csv(hand_path) if hand_path.exists() else pd.DataFrame()
    return runs, hand


def _infer_subtype_from_filename(filename: str) -> str:
    """VLP17_*.dm4 → VLP17, VLP20_*.dm4 → VLP20, VLP_100_*.dm3 → VLP_100, etc."""
    stem = Path(filename).stem
    # Match the leading non-digit prefix + first digit block
    import re
    m = re.match(r"^([A-Za-z_]+\d+)", stem)
    return m.group(1) if m else stem


def _flag_metric(value: float, ref: pd.Series, name: str,
                 lower_is_better: bool = False) -> Optional[dict]:
    """
    Return a flag dict if `value` is outside (median ± IQR_TOL × IQR) of `ref`,
    else None. Skips flagging when ref has fewer than MIN_N_FOR_FLAG points
    (returns an "info" dict instead so the user sees the comparison anyway).
    """
    ref = pd.Series(ref).dropna()
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if len(ref) == 0:
        return None
    median = float(ref.median())
    iqr    = float(ref.quantile(0.75) - ref.quantile(0.25))
    lo     = median - IQR_TOL * iqr
    hi     = median + IQR_TOL * iqr
    in_range = lo <= value <= hi
    severity = "info" if len(ref) < MIN_N_FOR_FLAG else ("ok" if in_range else "warn")
    if severity == "ok":
        return None
    delta = value - median
    return {
        "metric":    name,
        "value":     value,
        "ref_median": median,
        "ref_iqr":   iqr,
        "ref_n":     int(len(ref)),
        "delta":     delta,
        "severity":  severity,  # "warn" or "info"
    }


def _coerce_run_input(run: Any) -> dict:
    """Accept a `run()` dict, a Path to a run dir, or a Path to SUMMARY.md."""
    if isinstance(run, dict):
        return run
    p = Path(run)
    # Path to a directory containing SUMMARY.md / vlp_measurements.csv
    if p.is_dir():
        summary = p / "SUMMARY.md"
        csv     = p / "vlp_measurements.csv"
        if not csv.exists():
            raise FileNotFoundError(f"no vlp_measurements.csv under {p}")
        df = pd.read_csv(csv)
        # Reconstruct a minimal dict from the CSV — derives the same per-image
        # rows the in-process `run()` would have produced.
        from analysis.vlp_measure_v2 import _per_image_summary, _overall_summary, _script_version_with_tunables  # noqa: WPS433
        return {
            "sample_type":    "VLP",
            "batch_id":       p.name,
            "script_version": _script_version_with_tunables(),
            "summary":        _overall_summary(df),
            "per_image":      _per_image_summary(df),
            "outputs":        {"run_dir": str(p), "summary_md": str(summary)},
        }
    raise TypeError(f"unsupported run input: {type(run).__name__}")


# ── Layer 1: per-image intrinsic quality ──────────────────────────────────

# Metrics to check, with whether higher is "better" (purely cosmetic — eval
# is two-sided IQR, but the label helps in the report).
_LAYER1_METRICS = [
    ("wall_fit_success_rate", "higher better"),
    ("reliable_rate",         "higher better"),
    ("drop_rate",             "lower better"),
    ("median_wall_cv",        "lower better"),
    ("capsid_median_nm",      "neutral"),  # subtype-relative
    ("capsid_mean_nm",        "neutral"),
]


def _per_image_quality_check(per_image: list[dict], reference_runs: pd.DataFrame) -> list[dict]:
    """For each image in the new run, compare metrics against same-subtype reference."""
    rows: list[dict] = []
    for img in per_image:
        subtype = _infer_subtype_from_filename(img["filename"])
        ref_sub = reference_runs[reference_runs["sample_subtype"] == subtype]
        flags = []
        for metric, _direction in _LAYER1_METRICS:
            if metric not in img or metric not in ref_sub.columns:
                continue
            f = _flag_metric(img[metric], ref_sub[metric], metric)
            if f is not None:
                flags.append(f)
        rows.append({
            "filename":      img["filename"],
            "sample_subtype": subtype,
            "n_ref_images":  int(len(ref_sub)),
            "flags":         flags,
            "n_warns":       sum(1 for f in flags if f["severity"] == "warn"),
        })
    return rows


# ── Layer 2: sample-level hand-vs-script bias ─────────────────────────────

def _hand_aggregate(hand: pd.DataFrame, batch_id: str,
                    sample_subtype: Optional[str] = None) -> Optional[dict]:
    """
    Find hand data for this batch (preferred) or this sample_subtype (fallback).
    Returns {n, gold_mean_nm, capsid_mean_nm, source} or None.
    """
    if hand.empty:
        return None
    sub = hand[hand["batch_id"] == batch_id]
    source = "exact batch_id match"
    if sub.empty and sample_subtype is not None:
        # Fallback: hand entries whose batch_id contains the subtype.
        # (e.g. batch_id "VLP17_v2_initial" matches sample_subtype "VLP17")
        mask = hand["batch_id"].astype(str).str.contains(sample_subtype, regex=False)
        sub = hand[mask]
        source = f"sample_subtype fallback ({sample_subtype})"
    if sub.empty:
        return None
    cap = sub["hand_capsid_diameter_nm"].dropna()
    gold = sub["hand_gold_diameter_nm"].dropna() if "hand_gold_diameter_nm" in sub.columns else pd.Series(dtype=float)
    return {
        "n":              int(len(sub)),
        "n_capsid":       int(len(cap)),
        "n_gold":         int(len(gold)),
        "capsid_mean_nm": float(cap.mean()) if len(cap) else float("nan"),
        "capsid_std_nm":  float(cap.std())  if len(cap) >= 2 else float("nan"),
        "gold_mean_nm":   float(gold.mean()) if len(gold) else float("nan"),
        "gold_std_nm":    float(gold.std())  if len(gold) >= 2 else float("nan"),
        "source":         source,
    }


def _hand_vs_script_check(run_summary: dict, batch_id: str,
                          per_image: list[dict],
                          reference_hand: pd.DataFrame) -> Optional[dict]:
    """
    Compute hand−script delta for this batch (if hand data exists).
    Returns None if no hand data is available at all.
    """
    # Pick a sample_subtype representative for fallback. Use the most common
    # subtype across the new run's images.
    subtypes = [_infer_subtype_from_filename(p["filename"]) for p in per_image]
    common = max(set(subtypes), key=subtypes.count) if subtypes else None
    hand_agg = _hand_aggregate(reference_hand, batch_id, common)
    if hand_agg is None:
        return None
    script_capsid_mean = run_summary["summary"].get("capsid_mean_nm")
    script_gold_mean   = run_summary["summary"].get("gold_mean_nm")
    out = {
        "hand_n":              hand_agg["n"],
        "hand_source":         hand_agg["source"],
        "hand_capsid_mean_nm": hand_agg["capsid_mean_nm"],
        "hand_capsid_std_nm":  hand_agg["capsid_std_nm"],
        "script_capsid_mean_nm": script_capsid_mean,
        "delta_capsid_nm":     (hand_agg["capsid_mean_nm"] - script_capsid_mean)
                                  if script_capsid_mean is not None else float("nan"),
    }
    if not np.isnan(hand_agg["gold_mean_nm"]):
        out.update({
            "hand_gold_mean_nm":   hand_agg["gold_mean_nm"],
            "script_gold_mean_nm": script_gold_mean,
            "delta_gold_nm":       (hand_agg["gold_mean_nm"] - script_gold_mean)
                                      if script_gold_mean is not None else float("nan"),
        })
    return out


# ── Report rendering ──────────────────────────────────────────────────────

def _render_report(result: dict) -> str:
    s = result
    md = [
        f"# Eval — {s['batch_id']}",
        "",
        f"- **Sample:** {s['sample_type']}",
        f"- **Script:** `{s['script_version']}`",
        f"- **Reference:** `benchmarks/{s['sample_type'].lower()}/reference_runs.csv` "
        f"({s['n_ref_runs']} prior images)",
        "",
        "## Headline",
        f"- {s['n_images']} images evaluated",
        f"- **{s['n_warns_total']} per-image quality warning(s)**",
        f"- Hand vs script: " + (
            "no hand measurements registered for this batch" if s["hand_vs_script"] is None
            else f"Δ capsid {s['hand_vs_script']['delta_capsid_nm']:+.2f} nm "
                 f"(hand {s['hand_vs_script']['hand_capsid_mean_nm']:.2f} nm, "
                 f"script {s['hand_vs_script']['script_capsid_mean_nm']:.2f} nm; "
                 f"hand source: {s['hand_vs_script']['hand_source']})"
        ),
    ]

    md += ["", "## Per-image quality check",
           "", "*Compares each image's quality metrics against the same-subtype "
               "distribution in `reference_runs.csv`. Flags values outside median ± 2× IQR.*", ""]
    if any(r["flags"] for r in s["per_image_quality"]):
        md += [
            "| filename | subtype | n_ref | warnings |",
            "|---|---|---|---|",
        ]
        for r in s["per_image_quality"]:
            warns = "; ".join(
                f"`{f['metric']}` = {f['value']:.3f} "
                f"(ref median {f['ref_median']:.3f} ± IQR {f['ref_iqr']:.3f}; "
                f"Δ {f['delta']:+.3f})"
                for f in r["flags"] if f["severity"] == "warn"
            ) or "—"
            md.append(f"| {r['filename']} | {r['sample_subtype']} | {r['n_ref_images']} | {warns} |")
    else:
        md.append("All images within reference range. ✓")

    md += ["", "## Hand vs script",
           "", "*Compares this batch's aggregate measurements to hand-measured ground "
               "truth in `reference_hand.csv`.*", ""]
    if s["hand_vs_script"] is None:
        md.append("No hand measurements registered for this batch (or sample_subtype). "
                  "Drop a hand CSV into `results/<batch>/hand/` and re-run eval to enable.")
    else:
        H = s["hand_vs_script"]
        md.append(f"Hand source: **{H['hand_source']}**, n={H['hand_n']} particles.")
        md.append("")
        md.append("| metric | hand | script | Δ (hand − script) |")
        md.append("|---|---|---|---|")
        if "hand_gold_mean_nm" in H and not np.isnan(H["hand_gold_mean_nm"]):
            md.append(f"| Gold mean | {H['hand_gold_mean_nm']:.2f} nm "
                      f"| {H['script_gold_mean_nm']:.2f} nm "
                      f"| **{H['delta_gold_nm']:+.2f} nm** |")
        md.append(f"| Capsid mean | {H['hand_capsid_mean_nm']:.2f} nm "
                  f"| {H['script_capsid_mean_nm']:.2f} nm "
                  f"| **{H['delta_capsid_nm']:+.2f} nm** |")

    return "\n".join(md) + "\n"


# ── Public API ────────────────────────────────────────────────────────────

def evaluate(
    run: Any,
    *,
    sample_type: str = "VLP",
    benchmarks_dir: Path | str = DEFAULT_BENCH,
    write_report: bool = True,
) -> dict:
    """
    Compare a `run()` output (dict, run dir, or path to SUMMARY.md) against
    `benchmarks/<sample_type>/reference_runs.csv` and `reference_hand.csv`.

    Returns a dict with structured eval results and writes an `eval_report.md`
    next to the run's outputs (unless write_report=False).
    """
    benchmarks_dir = Path(benchmarks_dir)
    run_data       = _coerce_run_input(run)
    runs_ref, hand_ref = _load_reference(sample_type, benchmarks_dir)

    per_image_quality = _per_image_quality_check(run_data["per_image"], runs_ref)
    hand_vs_script    = _hand_vs_script_check(run_data, run_data["batch_id"], run_data["per_image"], hand_ref)

    n_warns_total = sum(r["n_warns"] for r in per_image_quality)

    result = {
        "sample_type":        sample_type,
        "batch_id":           run_data["batch_id"],
        "script_version":     run_data["script_version"],
        "n_images":           len(run_data["per_image"]),
        "n_ref_runs":         int(len(runs_ref)),
        "n_warns_total":      n_warns_total,
        "per_image_quality":  per_image_quality,
        "hand_vs_script":     hand_vs_script,
    }

    if write_report and "outputs" in run_data and "run_dir" in run_data["outputs"]:
        report_path = Path(run_data["outputs"]["run_dir"]) / "eval_report.md"
        report_path.write_text(_render_report(result))
        result["report_path"] = str(report_path)

    return result


# ── CLI ───────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run_dir", type=Path, help="Directory with vlp_measurements.csv")
    p.add_argument("--sample-type", default="VLP")
    p.add_argument("--benchmarks-dir", type=Path, default=DEFAULT_BENCH)
    p.add_argument("--no-report", action="store_true",
                   help="Don't write eval_report.md (still returns dict)")
    p.add_argument("--json", action="store_true",
                   help="Print full result as JSON instead of headline")
    args = p.parse_args()

    result = evaluate(args.run_dir, sample_type=args.sample_type,
                      benchmarks_dir=args.benchmarks_dir,
                      write_report=not args.no_report)
    if args.json:
        print(json.dumps(result, indent=2, default=str))
        return

    print(f"Eval — {result['batch_id']}")
    print(f"  reference:        {result['n_ref_runs']} prior runs")
    print(f"  per-image quality: {result['n_warns_total']} warning(s) across {result['n_images']} image(s)")
    if result["hand_vs_script"] is None:
        print("  hand vs script:   no hand measurements registered for this batch")
    else:
        H = result["hand_vs_script"]
        print(f"  hand vs script:   Δ capsid {H['delta_capsid_nm']:+.2f} nm "
              f"({H['hand_n']} hand particles, source: {H['hand_source']})")
    if "report_path" in result:
        print(f"  report:    {result['report_path']}")


if __name__ == "__main__":
    main()
