"""
eval.py — compare a measurement run against the reference benchmarks.

Two questions the eval answers:

  1. Per-image quality check
     Was each image measured WELL — clean fits, low noise, consistent across
     particles? This is asked using **size-invariant quality metrics only**
     (rates and CVs), not absolute particle sizes. Sizes are what the script
     measures; quality is what the eval judges.

     Metrics gated:
       wall_fit_success_rate, reliable_rate, drop_rate,
       median_wall_cv, iqr_wall_cv.

     Each is compared against the distribution of the SAME metric across all
     prior trusted runs in `benchmarks/<sample_type>/reference_runs.csv`
     (pooled — no subtype filter, because sample names are intent rather than
     measured size, and we don't want to lock the eval to mislabelled bins).

     Flag if the new value sits outside median ± 2 × IQR of the reference.

  2. Hand vs script
     If the scientist dropped one or more hand-measurement CSVs into
     `<run_dir>/hand/*.csv`, aggregate them and compute hand_mean − script_mean
     for capsid and gold. Surfaces calibration drift between the script and
     human ImageJ measurements.

     Filesystem-paired only — no batch_id, no subtype matching. If the
     scientist wants per-run hand validation, they put the CSV in the run's
     `hand/` folder.

Outputs:
  - structured dict (returned)
  - `<run_dir>/eval_report.md` (human-readable, with explanatory intro)

Usage:
  from analysis.eval import evaluate
  ev = evaluate(run_result_dict_or_run_dir, sample_type="VLP")

  # CLI:
  uv run python -m analysis.eval results/<sample_name>
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
MIN_N_FOR_FLAG  = 4      # need this many ref points before we'll flag (else just info)


# ── Reference loading ────────────────────────────────────────────────────

def _load_reference_runs(sample_type: str, benchmarks_dir: Path) -> pd.DataFrame:
    runs_path = benchmarks_dir / sample_type.lower() / "reference_runs.csv"
    if not runs_path.exists():
        raise FileNotFoundError(f"reference_runs.csv not found: {runs_path}")
    return pd.read_csv(runs_path)


# ── Per-image quality check ─────────────────────────────────────────────

# Size-invariant quality metrics. These describe HOW WELL we measured —
# regardless of absolute particle size. Adding metrics here means they get
# checked against the reference distribution; absolute-size metrics like
# capsid_median_nm or gold_mean_nm intentionally NOT here.
_QUALITY_METRICS = [
    "wall_fit_success_rate",
    "reliable_rate",
    "drop_rate",
    "median_wall_cv",
    "iqr_wall_cv",
]

_METRIC_DEFINITIONS = {
    "wall_fit_success_rate":
        "Fraction of detected gold particles for which the script found a clean "
        "capsid wall. Low = the protein ring isn't clear enough in this image "
        "(stain too thin, contrast too low, or focus issues).",
    "reliable_rate":
        "Fraction of detected gold particles whose capsid measurement passed "
        "all quality filters (good wall fit + circular wall). Low = many "
        "particles look problematic to the algorithm.",
    "drop_rate":
        "1 − reliable_rate. Fraction excluded from the final mean.",
    "median_wall_cv":
        "Per-particle, the script measures wall radius in 8 angular sectors "
        "and computes std/mean (= CV). 0 = perfect circle. The median CV "
        "across all particles in this image is reported. Higher = particles "
        "look non-circular (deformed walls, particle clustering, bad fits).",
    "iqr_wall_cv":
        "Spread of wall_cv across particles in this image. High IQR = some "
        "particles fit cleanly, others don't — often uneven staining.",
}


def _flag_metric(value: float, ref: pd.Series, name: str) -> Optional[dict]:
    """
    Return a flag dict if `value` is outside (median ± IQR_TOL × IQR) of `ref`,
    else None. Skips flagging when ref has fewer than MIN_N_FOR_FLAG points.
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
    return {
        "metric":     name,
        "value":      value,
        "ref_median": median,
        "ref_iqr":    iqr,
        "ref_n":      int(len(ref)),
        "delta":      value - median,
        "severity":   severity,  # "warn" or "info"
    }


def _per_image_quality_check(per_image: list[dict], reference_runs: pd.DataFrame) -> list[dict]:
    """For each image, compare each size-invariant metric against the pooled reference."""
    rows: list[dict] = []
    for img in per_image:
        flags = []
        for metric in _QUALITY_METRICS:
            if metric not in img or metric not in reference_runs.columns:
                continue
            f = _flag_metric(img[metric], reference_runs[metric], metric)
            if f is not None:
                flags.append(f)
        rows.append({
            "filename":     img["filename"],
            "n_ref_images": int(len(reference_runs)),
            "flags":        flags,
            "n_warns":      sum(1 for f in flags if f["severity"] == "warn"),
        })
    return rows


# ── Hand vs script (filesystem-paired) ───────────────────────────────────

def _read_hand_csvs(run_dir: Path) -> Optional[pd.DataFrame]:
    """
    Read any *.csv files under `<run_dir>/hand/`. Each CSV is the ImageJ
    output for some subset of particles in this run. Concatenate.

    Returns a DataFrame with at least `hand_capsid_diameter_nm` (per-particle),
    or None if no hand directory or no CSVs.

    The CSV format is whatever the scientist exported from ImageJ; we expect
    a `Length` column. Length unit (µm vs nm) is inferred from the magnitude:
    if max < 1.0 we treat as µm and convert; otherwise nm. Every existing
    hand CSV in this repo follows that convention.
    """
    hand_dir = run_dir / "hand"
    if not hand_dir.is_dir():
        return None
    csvs = sorted(hand_dir.glob("*.csv"))
    if not csvs:
        return None

    rows = []
    for csv_path in csvs:
        df = pd.read_csv(csv_path)
        if "Length" not in df.columns:
            continue
        lengths = df["Length"].to_numpy(dtype=float)
        if lengths.size == 0:
            continue
        unit = "um" if np.nanmax(lengths) < 1.0 else "nm"
        lengths_nm = lengths * 1000.0 if unit == "um" else lengths
        # Heuristic: even-row-count CSVs with Length variation across pairs are
        # most likely paired (gold, capsid) — same convention as seed data.
        # We treat odd-indexed rows (0-indexed: rows 1, 3, 5, …) as capsid.
        # Single-column hand CSVs (capsid only) get all rows used.
        if len(lengths_nm) % 2 == 0 and len(lengths_nm) >= 4:
            golds   = lengths_nm[0::2]
            capsids = lengths_nm[1::2]
            paired = True
        else:
            golds   = np.full_like(lengths_nm, np.nan)
            capsids = lengths_nm
            paired = False
        for i, (g, c) in enumerate(zip(golds, capsids), start=1):
            rows.append({
                "source_file":             csv_path.name,
                "particle_idx":            i,
                "hand_gold_diameter_nm":   float(g) if not np.isnan(g) else np.nan,
                "hand_capsid_diameter_nm": float(c),
                "paired_format":           paired,
            })
    if not rows:
        return None
    return pd.DataFrame(rows)


def _hand_vs_script_check(run_data: dict) -> Optional[dict]:
    """Read <run_dir>/hand/*.csv and compute hand−script delta. None if no CSVs."""
    run_dir_str = run_data.get("outputs", {}).get("run_dir")
    if not run_dir_str:
        return None
    hand = _read_hand_csvs(Path(run_dir_str))
    if hand is None or hand.empty:
        return None

    cap = hand["hand_capsid_diameter_nm"].dropna()
    gold = hand["hand_gold_diameter_nm"].dropna()
    n_csvs = hand["source_file"].nunique()

    script_capsid_mean = run_data["summary"].get("capsid_mean_nm")
    script_gold_mean   = run_data["summary"].get("gold_mean_nm")

    out = {
        "n_csvs":              int(n_csvs),
        "n_hand_particles":    int(len(hand)),
        "source_files":        sorted(hand["source_file"].unique().tolist()),
        "hand_capsid_mean_nm": float(cap.mean()) if len(cap) else float("nan"),
        "hand_capsid_std_nm":  float(cap.std())  if len(cap) >= 2 else float("nan"),
        "script_capsid_mean_nm": script_capsid_mean,
        "delta_capsid_nm":     (float(cap.mean()) - script_capsid_mean)
                                  if len(cap) and script_capsid_mean is not None else float("nan"),
    }
    if len(gold):
        out.update({
            "hand_gold_mean_nm":   float(gold.mean()),
            "hand_gold_std_nm":    float(gold.std()) if len(gold) >= 2 else float("nan"),
            "script_gold_mean_nm": script_gold_mean,
            "delta_gold_nm":       (float(gold.mean()) - script_gold_mean)
                                      if script_gold_mean is not None else float("nan"),
        })
    return out


# ── Run input coercion ───────────────────────────────────────────────────

def _coerce_run_input(run: Any) -> dict:
    """Accept a `run()` dict, a Path to a run dir, or a Path to SUMMARY.md."""
    if isinstance(run, dict):
        return run
    p = Path(run)
    if p.is_dir():
        csv = p / "vlp_measurements.csv"
        if not csv.exists():
            raise FileNotFoundError(f"no vlp_measurements.csv under {p}")
        df = pd.read_csv(csv)
        from analysis.vlp_measure_v2 import (  # noqa: WPS433
            _per_image_summary, _overall_summary, _script_version_with_tunables,
        )
        return {
            "sample_type":    "VLP",
            "sample_name":    p.name,
            "script_version": _script_version_with_tunables(),
            "summary":        _overall_summary(df),
            "per_image":      _per_image_summary(df),
            "outputs":        {"run_dir": str(p), "summary_md": str(p / "SUMMARY.md")},
        }
    raise TypeError(f"unsupported run input: {type(run).__name__}")


# ── Report rendering ─────────────────────────────────────────────────────

_INTRO = """\
> ## How to read this report
>
> This compares your run against the lab's reference set of past trusted runs.
> We check whether the **measurement quality** looks normal — *not* whether
> the absolute sizes are what you expected. Sizes are what the script
> measures; this report doesn't judge them. (Your VLP25 sample might turn out
> to be 20 nm. The eval doesn't care, as long as the measurement was clean.)
>
> ### What's a warning?
> Each per-image warning means a quality metric fell outside the typical range
> we see across past measurements. The "typical range" is **median ± 2 × IQR**
> (interquartile range — the spread of the middle 50% of past values).
> Roughly speaking, half the past runs sit within median ± 1 IQR; almost all
> within median ± 2 IQR. Outside that = unusual, worth a look.
>
> ### What to do when you see a warning
> Open the corresponding overlay PNG in `overlays/`. Most quality warnings
> show up visually — uneven staining, particle clustering, focus issues, weak
> contrast. The eye usually sees the cause within seconds. The numbers just
> tell you which image to look at.
>
> ### Metric reference
"""


def _render_report(result: dict) -> str:
    s = result
    md = [
        f"# Eval — {s['sample_name']}",
        "",
        f"- **Sample type:** {s['sample_type']}",
        f"- **Script:** `{s['script_version']}`",
        f"- **Reference:** `benchmarks/{s['sample_type'].lower()}/reference_runs.csv` "
        f"({s['n_ref_runs']} prior images)",
        "",
        _INTRO,
    ]
    for metric, definition in _METRIC_DEFINITIONS.items():
        md.append(f"> - **`{metric}`** — {definition}")
    md.append(">")
    md.append("> Absolute sizes (`capsid_mean_nm`, `gold_mean_nm`, etc.) are reported but "
              "**not gated**.")

    md += [
        "",
        "## Headline",
        f"- {s['n_images']} images evaluated",
        f"- **{s['n_warns_total']} per-image quality warning(s)**",
    ]
    sm = s.get("script_summary", {})
    if sm:
        md.append(
            f"- This run measured: gold {sm.get('gold_mean_nm', float('nan')):.2f} ± "
            f"{sm.get('gold_std_nm', float('nan')):.2f} nm, capsid "
            f"{sm.get('capsid_mean_nm', float('nan')):.2f} ± "
            f"{sm.get('capsid_std_nm', float('nan')):.2f} nm "
            f"(reported, not judged)"
        )
    if s["hand_vs_script"] is None:
        md.append("- Hand vs script: no hand CSVs in this run's `hand/` folder")
    else:
        H = s["hand_vs_script"]
        md.append(
            f"- Hand vs script: Δ capsid **{H['delta_capsid_nm']:+.2f} nm** "
            f"(hand n={H['n_hand_particles']} from {H['n_csvs']} CSV(s))"
        )

    md += ["", "## Per-image quality check", ""]
    if any(r["flags"] for r in s["per_image_quality"]):
        md += [
            "| filename | warnings |",
            "|---|---|",
        ]
        for r in s["per_image_quality"]:
            warns = "; ".join(
                f"`{f['metric']}` = {f['value']:.3f} "
                f"(ref median {f['ref_median']:.3f}, IQR {f['ref_iqr']:.3f}; "
                f"Δ {f['delta']:+.3f})"
                for f in r["flags"] if f["severity"] == "warn"
            )
            if warns:
                md.append(f"| `{r['filename']}` | {warns} |")
        md.append("")
        md.append(f"_Reference distribution built from {s['per_image_quality'][0]['n_ref_images']} "
                  "prior images of this sample type._")
    else:
        md.append("All images within reference range. ✓")

    md += ["", "## Hand vs script", ""]
    if s["hand_vs_script"] is None:
        md.append("No hand-measurement CSVs found in `<run_dir>/hand/`. Drop one or more "
                  "ImageJ CSVs there and re-run eval to see the script-vs-hand comparison.")
    else:
        H = s["hand_vs_script"]
        md.append(f"Read **{H['n_hand_particles']} particles** from "
                  f"{H['n_csvs']} CSV(s): {', '.join(H['source_files'])}")
        md.append("")
        md.append("| metric | hand | script | Δ (hand − script) |")
        md.append("|---|---|---|---|")
        if "hand_gold_mean_nm" in H:
            md.append(f"| Gold mean | {H['hand_gold_mean_nm']:.2f} nm "
                      f"| {H['script_gold_mean_nm']:.2f} nm "
                      f"| **{H['delta_gold_nm']:+.2f} nm** |")
        md.append(f"| Capsid mean | {H['hand_capsid_mean_nm']:.2f} nm "
                  f"| {H['script_capsid_mean_nm']:.2f} nm "
                  f"| **{H['delta_capsid_nm']:+.2f} nm** |")
        md.append("")
        md.append("_A persistent positive or negative delta across runs would suggest the "
                  "script reads systematically larger or smaller than human ImageJ tracing — "
                  "a calibration issue rather than a per-batch problem._")

    return "\n".join(md) + "\n"


# ── Public API ───────────────────────────────────────────────────────────

def evaluate(
    run: Any,
    *,
    sample_type: str = "VLP",
    benchmarks_dir: Path | str = DEFAULT_BENCH,
    write_report: bool = True,
) -> dict:
    """
    Compare a `run()` output (dict, run dir, or path to SUMMARY.md) against
    `benchmarks/<sample_type>/reference_runs.csv` plus any hand CSVs in
    `<run_dir>/hand/`.

    Returns a dict with structured eval results and writes an `eval_report.md`
    next to the run's outputs (unless write_report=False).
    """
    benchmarks_dir = Path(benchmarks_dir)
    run_data       = _coerce_run_input(run)
    runs_ref       = _load_reference_runs(sample_type, benchmarks_dir)

    per_image_quality = _per_image_quality_check(run_data["per_image"], runs_ref)
    hand_vs_script    = _hand_vs_script_check(run_data)
    n_warns_total     = sum(r["n_warns"] for r in per_image_quality)

    result = {
        "sample_type":        sample_type,
        "sample_name":        run_data.get("sample_name") or run_data.get("batch_id") or "(unnamed)",
        "script_version":     run_data["script_version"],
        "script_summary":     run_data.get("summary", {}),
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


# ── CLI ──────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("run_dir", type=Path)
    p.add_argument("--sample-type", default="VLP")
    p.add_argument("--benchmarks-dir", type=Path, default=DEFAULT_BENCH)
    p.add_argument("--no-report", action="store_true")
    p.add_argument("--json", action="store_true")
    args = p.parse_args()

    result = evaluate(args.run_dir, sample_type=args.sample_type,
                      benchmarks_dir=args.benchmarks_dir,
                      write_report=not args.no_report)
    if args.json:
        print(json.dumps(result, indent=2, default=str))
        return

    print(f"Eval — {result['sample_name']}")
    print(f"  reference:         {result['n_ref_runs']} prior runs")
    print(f"  per-image quality: {result['n_warns_total']} warning(s) across {result['n_images']} image(s)")
    if result["hand_vs_script"] is None:
        print("  hand vs script:    no hand CSVs in <run_dir>/hand/")
    else:
        H = result["hand_vs_script"]
        print(f"  hand vs script:    Δ capsid {H['delta_capsid_nm']:+.2f} nm "
              f"({H['n_hand_particles']} hand particles from {H['n_csvs']} CSV(s))")
    if "report_path" in result:
        print(f"  report:            {result['report_path']}")


if __name__ == "__main__":
    main()
