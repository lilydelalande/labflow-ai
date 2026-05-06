"""
add_to_reference.py — append approved runs / hand measurements to the
benchmark reference set. The manual gate that grows the eval baseline.

Two operations:

  add_run     — append per-image rows from a measurement run into
                `benchmarks/<sample_type>/reference_runs.csv`.

  add_hand    — append per-particle hand-measurement rows into
                `benchmarks/<sample_type>/reference_hand.csv`.

Both APPEND, never overwrite. The seed (`seed_benchmarks.py`) overwrites; this
script is for incremental growth after the initial seed.

Usage (CLI):
  uv run python -m analysis.add_to_reference run results/<batch_dir> \\
      --approver Lily --notes "good batch, hand-validated"

  uv run python -m analysis.add_to_reference hand path/to/hand.csv \\
      --batch-id VLP17_2026-05-15 --format paired --unit um \\
      --scientist Lily --measure-date 2026-05-15

Usage (programmatic):
  from analysis.add_to_reference import add_run, add_hand
  add_run(run_dir="results/foo", sample_type="VLP", approver="Lily", notes="...")
  add_hand(hand_csv="path/to/h.csv", batch_id="...", sample_type="VLP",
           hand_format="paired", length_unit="um", scientist="Lily")
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.seed_benchmarks import (
    _per_image_summary,
    _parse_paired_hand,
    _parse_single_hand,
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
            # force=True: drop matching rows from existing
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


def _infer_subtype_from_filename(filename: str) -> str:
    """VLP17_*.dm4 → VLP17, VLP_100_*.dm3 → VLP_100, etc."""
    import re
    stem = Path(filename).stem
    m = re.match(r"^([A-Za-z_]+\d+)", stem)
    return m.group(1) if m else stem


# ── add_run ───────────────────────────────────────────────────────────────

def add_run(
    run_dir: str | Path,
    *,
    sample_type: str = "VLP",
    batch_id: str | None = None,
    approver: str = "",
    notes: str = "",
    sample_subtype: str | None = None,
    benchmarks_dir: str | Path = DEFAULT_BENCH,
    force: bool = False,
) -> dict:
    """
    Append per-image rows from `run_dir/vlp_measurements.csv` to
    `benchmarks/<sample_type>/reference_runs.csv`.

    Dedup: refuses to add rows with the same (batch_id, filename) as existing
    entries unless `force=True`, in which case the existing rows are replaced.
    """
    run_dir = Path(run_dir)
    csv     = run_dir / "vlp_measurements.csv"
    if not csv.exists():
        raise FileNotFoundError(f"no vlp_measurements.csv under {run_dir}")
    if not approver:
        raise ValueError("approver is required (e.g. 'Lily')")

    bid = batch_id or run_dir.name

    rows = _per_image_summary(
        csv_path       = csv,
        batch_id       = bid,
        sample_subtype = sample_subtype or "MIXED",
        approver       = approver,
        notes          = notes,
    )
    if sample_subtype is None:
        rows["sample_subtype"] = rows["filename"].apply(_infer_subtype_from_filename)

    runs_csv, _ = _bench_paths(sample_type, Path(benchmarks_dir))
    n_before, n_after, n_replaced = _append_csv_dedup(
        rows, runs_csv,
        dedup_keys=["batch_id", "filename"],
        force=force,
    )

    return {
        "csv_path":   str(runs_csv),
        "n_added":    int(len(rows)),
        "n_replaced": n_replaced,
        "n_before":   n_before,
        "n_after":    n_after,
        "batch_id":   bid,
        "subtypes":   sorted(rows["sample_subtype"].unique().tolist()),
    }


# ── add_hand ──────────────────────────────────────────────────────────────

def add_hand(
    hand_csv: str | Path,
    *,
    batch_id: str,
    sample_type: str = "VLP",
    hand_format: str = "paired",     # "paired" (alternating gold/capsid) or "capsid_only"
    length_unit: str = "um",         # "um" or "nm"
    scientist: str = "",
    measure_date: str | None = None,
    image_filename: str | None = None,
    notes: str = "",
    benchmarks_dir: str | Path = DEFAULT_BENCH,
    force: bool = False,
) -> dict:
    """
    Append per-particle rows from a hand-measurement CSV to
    `benchmarks/<sample_type>/reference_hand.csv`.

    `hand_format`:
      - "paired"      — alternating gold (odd) / capsid (even) rows, 2N rows for N particles
      - "capsid_only" — every row is one capsid measurement (used for per-image diagnostics)

    `image_filename` is set on every appended row. Leave None for whole-batch
    aggregates; set to the .dm4 filename for per-image diagnostic measurements.
    """
    hand_csv = Path(hand_csv)
    if not hand_csv.exists():
        raise FileNotFoundError(f"hand CSV not found: {hand_csv}")
    if not scientist:
        raise ValueError("scientist is required (e.g. 'Lily')")

    measure_date = measure_date or _dt.date.today().isoformat()

    if hand_format == "paired":
        rows = _parse_paired_hand(
            csv_path     = hand_csv,
            length_unit  = length_unit,
            batch_id     = batch_id,
            scientist    = scientist,
            measure_date = measure_date,
            notes        = notes,
        )
    elif hand_format == "capsid_only":
        if image_filename is None:
            print("warning: capsid_only without image_filename — leaving image_filename=NULL",
                  file=sys.stderr)
        rows = _parse_single_hand(
            csv_path       = hand_csv,
            length_unit    = length_unit,
            batch_id       = batch_id,
            image_filename = image_filename or "",
            scientist      = scientist,
            measure_date   = measure_date,
            notes          = notes,
        )
        if image_filename is None:
            rows["image_filename"] = None
    else:
        raise ValueError(f"unknown hand_format: {hand_format!r} (expected 'paired' or 'capsid_only')")

    if image_filename is not None:
        rows["image_filename"] = image_filename

    _, hand_path = _bench_paths(sample_type, Path(benchmarks_dir))
    # Dedup on (batch_id, source_file): the same hand CSV uploaded twice for
    # the same batch is the duplicate case we're guarding against.
    n_before, n_after, n_replaced = _append_csv_dedup(
        rows, hand_path,
        dedup_keys=["batch_id", "source_file"],
        force=force,
    )

    return {
        "csv_path":       str(hand_path),
        "n_added":        int(len(rows)),
        "n_replaced":     n_replaced,
        "n_before":       n_before,
        "n_after":        n_after,
        "batch_id":       batch_id,
        "image_filename": image_filename,
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
    pr.add_argument("run_dir", type=Path,
                    help="Directory containing vlp_measurements.csv")
    pr.add_argument("--sample-type", default="VLP")
    pr.add_argument("--batch-id", default=None,
                    help="Defaults to the folder name of run_dir")
    pr.add_argument("--approver", required=True,
                    help="Person approving this run (recorded in CSV)")
    pr.add_argument("--notes", default="",
                    help="Free-form note (why this run is reference-worthy)")
    pr.add_argument("--sample-subtype", default=None,
                    help="Force a single subtype (VLP17, VLP20, ...). "
                         "Default: auto-infer per image from filename.")
    pr.add_argument("--force", action="store_true",
                    help="If (batch_id, filename) rows already exist in "
                         "reference_runs.csv, REPLACE them instead of refusing.")
    pr.add_argument("--benchmarks-dir", type=Path, default=DEFAULT_BENCH)

    ph = sub.add_parser("hand", help="Append hand measurements to reference_hand.csv")
    ph.add_argument("hand_csv", type=Path)
    ph.add_argument("--batch-id", required=True,
                    help="Join key. Use the same batch_id as the script run "
                         "this hand data validates.")
    ph.add_argument("--sample-type", default="VLP")
    ph.add_argument("--format", dest="hand_format",
                    choices=["paired", "capsid_only"], default="paired",
                    help="paired = alternating gold/capsid rows; "
                         "capsid_only = every row is a capsid measurement")
    ph.add_argument("--unit", dest="length_unit", choices=["um", "nm"], default="um",
                    help="Length unit in the CSV's 'Length' column")
    ph.add_argument("--scientist", required=True)
    ph.add_argument("--measure-date", default=None,
                    help="ISO date (default: today)")
    ph.add_argument("--image-filename", default=None,
                    help="Per-image diagnostic? Pass the .dm4 filename. "
                         "Whole-batch aggregate? Leave unset.")
    ph.add_argument("--notes", default="")
    ph.add_argument("--force", action="store_true",
                    help="If (batch_id, source_file) rows already exist in "
                         "reference_hand.csv, REPLACE them instead of refusing.")
    ph.add_argument("--benchmarks-dir", type=Path, default=DEFAULT_BENCH)

    args = p.parse_args()

    try:
        if args.cmd == "run":
            result = add_run(
                run_dir         = args.run_dir,
                sample_type     = args.sample_type,
                batch_id        = args.batch_id,
                approver        = args.approver,
                notes           = args.notes,
                sample_subtype  = args.sample_subtype,
                benchmarks_dir  = args.benchmarks_dir,
                force           = args.force,
            )
            _format_result("appended run rows", result)
        elif args.cmd == "hand":
            result = add_hand(
                hand_csv       = args.hand_csv,
                batch_id       = args.batch_id,
                sample_type    = args.sample_type,
                hand_format    = args.hand_format,
                length_unit    = args.length_unit,
                scientist      = args.scientist,
                measure_date   = args.measure_date,
                image_filename = args.image_filename,
                notes          = args.notes,
                benchmarks_dir = args.benchmarks_dir,
                force          = args.force,
            )
            _format_result("appended hand rows", result)
    except DuplicateReferenceError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
