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

import os

# Force headless matplotlib backend before any worker imports vlp_measure_v2
# / bmv_measure (which import pyplot). Avoids a Python/Tk dock icon flashing
# on macOS for every parallel worker, and avoids GUI-backend hangs in CI.
os.environ.setdefault("MPLBACKEND", "Agg")

# Pin numerical-library thread counts to 1 so worker subprocesses are
# *deterministic*. Without this, scikit-image's Hough circle transform
# and CLAHE multithread internally; with 4 worker processes × N threads
# each, accumulator-order non-determinism flips borderline particles'
# quality classification and gives different reliable_rate from run to
# run. Each worker stays single-threaded; multiple workers stay parallel.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
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


def _validate_one_row(args: tuple) -> dict:
    """Worker fn — re-measure one image and return its row dict.
    Top-level so ProcessPoolExecutor can pickle it."""
    sample_type, sample_name, filename, img_path_str, ref_dict, tolerances, bmv_expected_nm = args
    img_path = Path(img_path_str)
    t0 = time.monotonic()
    try:
        actual = _measure_one_image(img_path, sample_type=sample_type,
                                     expected_nm=bmv_expected_nm)
    except Exception as exc:
        return {"sample_name": sample_name, "filename": filename,
                "status": "error", "error": str(exc), "drifts": [],
                "elapsed_s": time.monotonic() - t0}
    elapsed = time.monotonic() - t0
    ref = pd.Series(ref_dict)
    drifts = _compare(ref, actual, tolerances)
    return {"sample_name": sample_name, "filename": filename,
            "status": ("ok" if not drifts else "drift"),
            "drifts": drifts, "actual": actual, "elapsed_s": elapsed}


def _print_row_progress(row: dict) -> None:
    """Stream a one-line progress note as each row finishes."""
    label = f"{row['sample_name']}/{row['filename']}"
    elapsed = row.get("elapsed_s")
    timing = f"  [{elapsed:5.1f}s]" if elapsed is not None else ""
    if row["status"] == "ok":
        print(f"  ✓{timing} {label}", flush=True)
    elif row["status"] == "drift":
        print(f"  ✗{timing} {label}  REGRESSION", flush=True)
    elif row["status"] == "error":
        print(f"  !{timing} {label}  ERROR: {row.get('error')}", flush=True)


def validate(
    sample_type: str,
    references_dir: Path,
    benchmarks_dir: Path = DEFAULT_BENCH,
    tolerances: dict[str, float] | None = None,
    bmv_expected_nm: float = 28.0,
    workers: int = 4,
) -> dict:
    """Re-measure every reference row of one sample_type and report drifts."""
    tolerances = tolerances or DEFAULT_TOLERANCES
    runs_path = benchmarks_dir / sample_type.lower() / "reference_runs.csv"
    if not runs_path.exists():
        return {"sample_type": sample_type, "status": "no-reference",
                "reason": f"{runs_path} not found"}

    refs = pd.read_csv(runs_path)
    rows = []
    work = []  # (sample_type, sample_name, filename, img_path, ref_dict, tolerances, bmv_expected_nm)
    for _, ref in refs.iterrows():
        sample_name = ref["sample_name"]
        filename    = ref["filename"]
        img_path    = references_dir / sample_type.lower() / sample_name / filename
        if not img_path.exists():
            rows.append({"sample_name": sample_name, "filename": filename,
                         "status": "missing-image", "image_path": str(img_path),
                         "drifts": []})
            continue
        work.append((sample_type, sample_name, filename, str(img_path),
                     ref.to_dict(), tolerances, bmv_expected_nm))

    print(f"  ── {sample_type}: re-measuring {len(work)} reference rows "
          f"({len(refs) - len(work)} skipped — no source image) "
          f"with {min(workers, max(len(work), 1))} parallel worker(s)…", flush=True)

    if workers > 1 and len(work) > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = [ex.submit(_validate_one_row, w) for w in work]
            for fut in as_completed(futures):
                row = fut.result()
                rows.append(row)
                _print_row_progress(row)
    else:
        for w in work:
            row = _validate_one_row(w)
            rows.append(row)
            _print_row_progress(row)

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
        n_measured = res["n_ok"] + res["n_drift"]
        print(f"── {res['sample_type']}: {n_measured} of {res['n_ref']} reference rows re-measured  "
              f"(✓ {res['n_ok']}  ✗ {res['n_drift']}  ! {res['n_error']})")
        # Per-row entries only for cases that need attention. Healthy rows
        # are summarised by the count above.
        for r in res["rows"]:
            label = f"{r['sample_name']}/{r['filename']}"
            if r["status"] == "ok":
                print(f"  ✓ {label}")
            elif r["status"] == "drift":
                any_drift = True
                detail = "; ".join(
                    f"{d['metric']} {d['expected']:.3f} → {d['actual']:.3f} "
                    f"(Δ {d['delta']:+.3f}, tol {d['tolerance']:.3f})"
                    for d in r["drifts"]
                )
                print(f"  ✗ {label}  REGRESSION: {detail}")
            elif r["status"] == "error":
                print(f"  ! {label}  ERROR: {r.get('error')}")
            # missing-image: deliberately silent at the row level — counted below
        if res["n_missing"]:
            print(f"        {res['n_missing']} row(s) have no source image in the references repo")
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
    p.add_argument("--workers", type=int, default=4,
                   help="Parallel worker processes (default 4). Set to 1 for serial.")

    # Per-metric tolerances
    p.add_argument("--tol-capsid-nm",      type=float, default=DEFAULT_TOLERANCES["capsid_median_nm"])
    p.add_argument("--tol-wall-fit-rate",  type=float, default=DEFAULT_TOLERANCES["wall_fit_success_rate"])
    p.add_argument("--tol-reliable-rate",  type=float, default=DEFAULT_TOLERANCES["reliable_rate"])
    p.add_argument("--tol-wall-cv",        type=float, default=DEFAULT_TOLERANCES["median_wall_cv"])

    args = p.parse_args()
    print(f"validate_script starting — references: {args.references_dir} | workers: {args.workers}", flush=True)
    tolerances = {
        "capsid_median_nm":      args.tol_capsid_nm,
        "wall_fit_success_rate": args.tol_wall_fit_rate,
        "reliable_rate":         args.tol_reliable_rate,
        "median_wall_cv":        args.tol_wall_cv,
    }

    sample_types = ["VLP", "BMV"] if args.sample_type == "all" else [args.sample_type]
    t_start = time.monotonic()
    results = [
        validate(st, references_dir=args.references_dir,
                 benchmarks_dir=args.benchmarks_dir,
                 tolerances=tolerances, bmv_expected_nm=args.bmv_expected_nm,
                 workers=args.workers)
        for st in sample_types
    ]
    total_elapsed = time.monotonic() - t_start
    code = _print_report(results)

    print(f"Tolerances: {tolerances}")
    print(f"Total wall-clock: {total_elapsed:.1f}s")
    if code != 0:
        print("\nFAILED — regressions found above.", file=sys.stderr)
    else:
        print("\nClean — no regressions.")
    sys.exit(code)


if __name__ == "__main__":
    main()
