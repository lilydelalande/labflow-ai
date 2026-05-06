"""
add_to_reference.py — append approved runs / hand measurements to the
benchmark reference set. The manual gate that grows the eval baseline.

Two operations:

  add_run     — append per-image rows from a measurement run into
                `benchmarks/<sample_type>/reference_runs.csv`.
                Dedup keys: (sample_name, filename).

  add_hand    — append per-particle hand-measurement rows into
                `benchmarks/<sample_type>/reference_hand.csv`. Optionally also
                copies the CSV into `<run_dir>/hand/` so the per-batch eval
                script-vs-hand comparison fires on subsequent eval runs.
                Dedup keys: (sample_name, source_file).

Both APPEND, never overwrite. Use `force=True` / `--force` to replace
existing rows that collide on the dedup keys.

Usage (CLI):
  uv run python -m analysis.add_to_reference run results/<sample_name> \\
      --approver Lily --notes "good batch, hand-validated"

  uv run python -m analysis.add_to_reference hand path/to/hand.csv \\
      --sample-name VLP17_2026-05-15 --unit um \\
      --scientist Lily --run-dir results/<sample_name>

Usage (programmatic):
  from analysis.add_to_reference import add_run, add_hand
  add_run(run_dir="results/foo", sample_type="VLP", approver="Lily", notes="…")
  add_hand(hand_csv="path/h.csv", sample_name="…", sample_type="VLP",
           length_unit="um", scientist="Lily", run_dir="results/foo")
"""

from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.seed_benchmarks import (
    _per_image_summary,
    _parse_paired_hand,
    _safe_relpath,
)

REPO_ROOT     = Path(__file__).resolve().parent.parent
DEFAULT_BENCH = REPO_ROOT / "benchmarks"


# ── Helpers ───────────────────────────────────────────────────────────────

def _bench_paths(sample_type: str, benchmarks_dir: Path) -> tuple[Path, Path]:
    sub = benchmarks_dir / sample_type.lower()
    sub.mkdir(parents=True, exist_ok=True)
    return sub / "reference_runs.csv", sub / "reference_hand.csv"


class DuplicateReferenceError(ValueError):
    """Raised when add_run / add_hand would create duplicate rows. Use force=True to override."""


def _append_csv_dedup(
    new_rows: pd.DataFrame,
    csv_path: Path,
    *,
    dedup_keys: list[str],
    force: bool,
) -> tuple[int, int, int]:
    """
    Append `new_rows` to `csv_path`. Refuses if any row collides with an
    existing row on `dedup_keys`. With `force=True`, removes the colliding
    existing rows first (so the new rows replace them).

    Returns (n_before, n_after, n_replaced). Raises DuplicateReferenceError
    if duplicates exist and force=False.
    """
    if csv_path.exists():
        existing = pd.read_csv(csv_path)
    else:
        existing = pd.DataFrame()

    n_before = len(existing)

    if not existing.empty and all(k in existing.columns for k in dedup_keys):
        keys_existing = set(map(tuple, existing[dedup_keys].itertuples(index=False, name=None)))
        keys_new      = set(map(tuple, new_rows[dedup_keys].itertuples(index=False, name=None)))
        overlap = keys_existing & keys_new
        if overlap:
            sample = list(overlap)[:5]
            sample_str = "\n    ".join(
                ", ".join(f"{k}={v!r}" for k, v in zip(dedup_keys, t))
                for t in sample
            )
            if not force:
                raise DuplicateReferenceError(
                    f"refusing to append {len(overlap)} row(s) that already exist in "
                    f"{csv_path.name} (dedup keys: {dedup_keys}).\n"
                    f"  examples:\n    {sample_str}\n"
                    f"  pass force=True (or --force on the CLI) to REPLACE the existing rows."
                )
            mask = pd.Series(False, index=existing.index)
            for t in overlap:
                cond = pd.Series(True, index=existing.index)
                for k, v in zip(dedup_keys, t):
                    cond &= (existing[k] == v)
                mask |= cond
            n_replaced = int(mask.sum())
            existing   = existing[~mask].reset_index(drop=True)
        else:
            n_replaced = 0
    else:
        n_replaced = 0

    merged = pd.concat([existing, new_rows], ignore_index=True) if not existing.empty else new_rows
    merged.to_csv(csv_path, index=False)
    return n_before, len(merged), n_replaced


# ── add_run ───────────────────────────────────────────────────────────────

def add_run(
    run_dir: str | Path,
    *,
    sample_type: str = "VLP",
    sample_name: str | None = None,
    approver: str = "",
    notes: str = "",
    benchmarks_dir: str | Path = DEFAULT_BENCH,
    force: bool = False,
) -> dict:
    """
    Append per-image rows from `run_dir/vlp_measurements.csv` to
    `benchmarks/<sample_type>/reference_runs.csv`.

    `sample_name` defaults to the run_dir folder name. Dedup on
    (sample_name, filename); pass force=True to replace existing rows.
    """
    run_dir = Path(run_dir)
    csv     = run_dir / "vlp_measurements.csv"
    if not csv.exists():
        raise FileNotFoundError(f"no vlp_measurements.csv under {run_dir}")
    if not approver:
        raise ValueError("approver is required (e.g. 'Lily')")

    name = sample_name or run_dir.name
    rows = _per_image_summary(
        csv_path    = csv,
        sample_name = name,
        approver    = approver,
        notes       = notes,
    )

    runs_csv, _ = _bench_paths(sample_type, Path(benchmarks_dir))
    n_before, n_after, n_replaced = _append_csv_dedup(
        rows, runs_csv,
        dedup_keys=["sample_name", "filename"],
        force=force,
    )
    return {
        "csv_path":   str(runs_csv),
        "n_added":    int(len(rows)),
        "n_replaced": n_replaced,
        "n_before":   n_before,
        "n_after":    n_after,
        "sample_name": name,
    }


# ── add_hand ──────────────────────────────────────────────────────────────

def add_hand(
    hand_csv: str | Path,
    *,
    sample_name: str,
    sample_type: str = "VLP",
    length_unit: str = "um",
    scientist: str = "",
    measure_date: str | None = None,
    notes: str = "",
    run_dir: str | Path | None = None,
    benchmarks_dir: str | Path = DEFAULT_BENCH,
    force: bool = False,
) -> dict:
    """
    Append per-particle rows from a paired (gold + capsid) ImageJ CSV to
    `benchmarks/<sample_type>/reference_hand.csv`.

    If `run_dir` is provided, ALSO copy the CSV into `<run_dir>/hand/` so
    `evaluate(<run_dir>)` picks it up for the per-run hand-vs-script section.

    Dedup on (sample_name, source_file); force=True to replace existing rows.
    """
    hand_csv = Path(hand_csv)
    if not hand_csv.exists():
        raise FileNotFoundError(f"hand CSV not found: {hand_csv}")
    if not scientist:
        raise ValueError("scientist is required (e.g. 'Lily')")

    measure_date = measure_date or _dt.date.today().isoformat()

    # Optionally also drop the CSV into <run_dir>/hand/ for filesystem-paired eval.
    copied_to: str | None = None
    if run_dir is not None:
        run_dir_path = Path(run_dir)
        hand_dir = run_dir_path / "hand"
        hand_dir.mkdir(parents=True, exist_ok=True)
        target = hand_dir / hand_csv.name
        if target.resolve() != hand_csv.resolve():
            shutil.copy2(hand_csv, target)
        copied_to = str(target)

    rows = _parse_paired_hand(
        csv_path     = hand_csv,
        length_unit  = length_unit,
        sample_name  = sample_name,
        scientist    = scientist,
        measure_date = measure_date,
        notes        = notes,
    )

    _, hand_path = _bench_paths(sample_type, Path(benchmarks_dir))
    n_before, n_after, n_replaced = _append_csv_dedup(
        rows, hand_path,
        dedup_keys=["sample_name", "source_file"],
        force=force,
    )
    return {
        "csv_path":       str(hand_path),
        "n_added":        int(len(rows)),
        "n_replaced":     n_replaced,
        "n_before":       n_before,
        "n_after":        n_after,
        "sample_name":    sample_name,
        "copied_to_run":  copied_to,
        "capsid_mean_nm": float(rows["hand_capsid_diameter_nm"].mean()),
        "gold_mean_nm":   float(rows["hand_gold_diameter_nm"].mean())
                              if rows["hand_gold_diameter_nm"].notna().any() else None,
    }


# ── CLI ───────────────────────────────────────────────────────────────────

def _format_result(label: str, r: dict) -> None:
    print(f"\n── {label} ──")
    for k, v in r.items():
        print(f"  {k}: {v}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pr = sub.add_parser("run", help="Append a measurement run to reference_runs.csv")
    pr.add_argument("run_dir", type=Path)
    pr.add_argument("--sample-type", default="VLP")
    pr.add_argument("--sample-name", default=None,
                    help="Defaults to the folder name of run_dir")
    pr.add_argument("--approver", required=True)
    pr.add_argument("--notes", default="")
    pr.add_argument("--force", action="store_true")
    pr.add_argument("--benchmarks-dir", type=Path, default=DEFAULT_BENCH)

    ph = sub.add_parser("hand", help="Append hand measurements to reference_hand.csv")
    ph.add_argument("hand_csv", type=Path)
    ph.add_argument("--sample-name", required=True,
                    help="Identifier shared with the script run this hand data validates "
                         "(usually the run folder name).")
    ph.add_argument("--sample-type", default="VLP")
    ph.add_argument("--unit", dest="length_unit", choices=["um", "nm"], default="um")
    ph.add_argument("--scientist", required=True)
    ph.add_argument("--measure-date", default=None)
    ph.add_argument("--notes", default="")
    ph.add_argument("--run-dir", type=Path, default=None,
                    help="Also copy the hand CSV into <run-dir>/hand/ so eval can pair it.")
    ph.add_argument("--force", action="store_true")
    ph.add_argument("--benchmarks-dir", type=Path, default=DEFAULT_BENCH)

    args = p.parse_args()

    try:
        if args.cmd == "run":
            result = add_run(
                run_dir         = args.run_dir,
                sample_type     = args.sample_type,
                sample_name     = args.sample_name,
                approver        = args.approver,
                notes           = args.notes,
                benchmarks_dir  = args.benchmarks_dir,
                force           = args.force,
            )
            _format_result("appended run rows", result)
        elif args.cmd == "hand":
            result = add_hand(
                hand_csv       = args.hand_csv,
                sample_name    = args.sample_name,
                sample_type    = args.sample_type,
                length_unit    = args.length_unit,
                scientist      = args.scientist,
                measure_date   = args.measure_date,
                notes          = args.notes,
                run_dir        = args.run_dir,
                benchmarks_dir = args.benchmarks_dir,
                force          = args.force,
            )
            _format_result("appended hand rows", result)
    except DuplicateReferenceError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
