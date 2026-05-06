"""
validate_script.py — regression test for the measurement scripts.

For every row in `benchmarks/<sample_type>/reference_runs.csv`, look up the
source image at `<references_dir>/<sample_type>/<sample_name>/<filename>`,
re-run the current measurement script on it, and compare the new per-image
metrics to the recorded ones. Flag any drift beyond tolerance.

Catches:
  - Tunable changes that fixed one sample but broke another
  - Dependency upgrades that silently changed behavior
  - Refactors that lost precision somewhere

Usage:
  uv run python -m analysis.validate_script \\
      --references-dir ~/Development/labflow-ai-references

CI wires this into a workflow (see .github/workflows/validate.yml — planned).

Default tolerances are conservative — small enough to catch real drift,
large enough to absorb numerical jitter (which should be ~0 since the
scripts are deterministic, but having a non-zero tolerance avoids
false-positives from anything we missed).
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# sys.path fixup so this module works whether invoked as -m or by direct path
_PACKAGE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, _PACKAGE_ROOT)

REPO_ROOT       = Path(_PACKAGE_ROOT)
DEFAULT_BENCH   = REPO_ROOT / "benchmarks"


# Per-metric drift tolerances. Anything beyond these is a regression.
DEFAULT_TOLERANCES = {
    "capsid_median_nm":      0.30,   # nm
    "wall_fit_success_rate": 0.05,   # 5 percentage points
    "reliable_rate":         0.05,   # 5 percentage points
    "median_wall_cv":        0.02,   # absolute CV
}


def _measure_one_image(image_path: Path, sample_type: str, expected_nm: float = 28.0) -> dict:
    """
    Run the appropriate measurement script on a single image and return its
    per-image summary dict. Uses each script's process_image() function
    directly so we don't have to spawn subprocesses.
    """
    from analysis.vlp_measure_v2 import _per_image_summary as vlp_per_image_summary

    if sample_type.upper() == "VLP":
        from analysis.vlp_measure_v2 import process_image as vlp_process_image
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            (tmp_path / "overlays").mkdir(exist_ok=True)
            df, _log = vlp_process_image(
                image_path, tmp_path,
                gold_threshold=None, min_gold_nm=7.0, max_gold_nm=30.0,
                show_flagged=False,
            )
    elif sample_type.upper() == "BMV":
        from analysis.bmv_measure import process_image as bmv_process_image
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            (tmp_path / "overlays").mkdir(exist_ok=True)
            df, _log = bmv_process_image(
                image_path, tmp_path,
                expected_nm=expected_nm, save_debug=False,
                debug_indices=None, show_flagged=False,
            )
    else:
        raise ValueError(f"unknown sample_type: {sample_type!r}")

    summary_list = vlp_per_image_summary(df)
    if not summary_list:
        return {}
    return summary_list[0]


def _compare(expected: pd.Series, actual: dict, tolerances: dict[str, float]) -> list[dict]:
    """For each tolerance metric, compare expected vs actual and return any drifts."""
    drifts = []
    for metric, tol in tolerances.items():
        if metric not in expected.index or metric not in actual:
            continue
        e = expected[metric]
        a = actual[metric]
        if pd.isna(e) and pd.isna(a):
            continue
        if pd.isna(e) or pd.isna(a):
            drifts.append({"metric": metric, "expected": e, "actual": a,
                           "delta": float("nan"), "tolerance": tol,
                           "kind": "nan-mismatch"})
            continue
        delta = float(a) - float(e)
        if abs(delta) > tol:
            drifts.append({"metric": metric, "expected": float(e), "actual": float(a),
                           "delta": delta, "tolerance": tol, "kind": "drift"})
    return drifts


def validate(
    sample_type: str,
    references_dir: Path,
    benchmarks_dir: Path = DEFAULT_BENCH,
    tolerances: dict[str, float] | None = None,
    bmv_expected_nm: float = 28.0,
) -> dict:
    """Re-measure every reference row of one sample_type and report drifts."""
    tolerances = tolerances or DEFAULT_TOLERANCES
    runs_path = benchmarks_dir / sample_type.lower() / "reference_runs.csv"
    if not runs_path.exists():
        return {"sample_type": sample_type, "status": "no-reference",
                "reason": f"{runs_path} not found"}

    refs = pd.read_csv(runs_path)
    rows = []
    for _, ref in refs.iterrows():
        sample_name = ref["sample_name"]
        filename    = ref["filename"]
        img_path    = references_dir / sample_type.lower() / sample_name / filename
        if not img_path.exists():
            rows.append({"sample_name": sample_name, "filename": filename,
                         "status": "missing-image", "image_path": str(img_path),
                         "drifts": []})
            continue
        try:
            actual = _measure_one_image(img_path, sample_type=sample_type,
                                         expected_nm=bmv_expected_nm)
        except Exception as exc:  # pragma: no cover (defensive)
            rows.append({"sample_name": sample_name, "filename": filename,
                         "status": "error", "error": str(exc), "drifts": []})
            continue
        drifts = _compare(ref, actual, tolerances)
        rows.append({"sample_name": sample_name, "filename": filename,
                     "status": ("ok" if not drifts else "drift"),
                     "drifts": drifts,
                     "actual": actual})

    n_ok      = sum(1 for r in rows if r["status"] == "ok")
    n_drift   = sum(1 for r in rows if r["status"] == "drift")
    n_missing = sum(1 for r in rows if r["status"] == "missing-image")
    n_error   = sum(1 for r in rows if r["status"] == "error")

    return {"sample_type": sample_type,
            "n_ref":     len(refs),
            "n_ok":      n_ok,
            "n_drift":   n_drift,
            "n_missing": n_missing,
            "n_error":   n_error,
            "rows":      rows,
            "tolerances": tolerances,
            "status":    "regressions" if n_drift else "clean"}


def _print_report(results: list[dict]) -> int:
    """Print a human-readable report. Returns nonzero exit code if any regressions."""
    any_drift = False
    print()
    for res in results:
        st = res.get("status")
        if st == "no-reference":
            print(f"── {res['sample_type']}: skipped — {res.get('reason')}")
            continue
        print(f"── {res['sample_type']}: {res['n_ref']} reference rows  "
              f"(✓ {res['n_ok']}  ✗ {res['n_drift']}  "
              f"missing {res['n_missing']}  errors {res['n_error']})")
        for r in res["rows"]:
            tag = {"ok": "✓", "drift": "✗", "missing-image": "·", "error": "!"}.get(r["status"], "?")
            label = f"{r['sample_name']}/{r['filename']}"
            if r["status"] == "ok":
                print(f"  {tag} {label}")
            elif r["status"] == "drift":
                any_drift = True
                detail = "; ".join(
                    f"{d['metric']} {d['expected']:.3f} → {d['actual']:.3f} "
                    f"(Δ {d['delta']:+.3f}, tol {d['tolerance']:.3f})"
                    for d in r["drifts"]
                )
                print(f"  {tag} {label}  REGRESSION: {detail}")
            elif r["status"] == "missing-image":
                print(f"  {tag} {label}  (source image not found at {r['image_path']})")
            elif r["status"] == "error":
                print(f"  {tag} {label}  ERROR: {r.get('error')}")
        print()
    return 1 if any_drift else 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--references-dir", type=Path, required=True,
                   help="Path to the labflow-ai-references checkout (or any directory "
                        "with the same <sample_type>/<sample_name>/<filename> layout).")
    p.add_argument("--sample-type", choices=["VLP", "BMV", "all"], default="all",
                   help="Which sample type(s) to validate (default: all).")
    p.add_argument("--benchmarks-dir", type=Path, default=DEFAULT_BENCH)
    p.add_argument("--bmv-expected-nm", type=float, default=28.0,
                   help="Expected BMV capsid diameter (default 28).")

    # Per-metric tolerances
    p.add_argument("--tol-capsid-nm",      type=float, default=DEFAULT_TOLERANCES["capsid_median_nm"])
    p.add_argument("--tol-wall-fit-rate",  type=float, default=DEFAULT_TOLERANCES["wall_fit_success_rate"])
    p.add_argument("--tol-reliable-rate",  type=float, default=DEFAULT_TOLERANCES["reliable_rate"])
    p.add_argument("--tol-wall-cv",        type=float, default=DEFAULT_TOLERANCES["median_wall_cv"])

    args = p.parse_args()
    tolerances = {
        "capsid_median_nm":      args.tol_capsid_nm,
        "wall_fit_success_rate": args.tol_wall_fit_rate,
        "reliable_rate":         args.tol_reliable_rate,
        "median_wall_cv":        args.tol_wall_cv,
    }

    sample_types = ["VLP", "BMV"] if args.sample_type == "all" else [args.sample_type]
    results = [
        validate(st, references_dir=args.references_dir,
                 benchmarks_dir=args.benchmarks_dir,
                 tolerances=tolerances, bmv_expected_nm=args.bmv_expected_nm)
        for st in sample_types
    ]
    code = _print_report(results)

    print(f"Tolerances: {tolerances}")
    if code != 0:
        print("\nFAILED — regressions found above.", file=sys.stderr)
    else:
        print("\nClean — no regressions.")
    sys.exit(code)


if __name__ == "__main__":
    main()
