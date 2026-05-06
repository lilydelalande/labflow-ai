"""
seed_benchmarks.py — Initial seed of benchmarks/<sample_type>/{reference_runs.csv,reference_hand.csv}.

Reads the existing v2 outputs in `results/vlp{17,20,100}_v2/` and the hand
measurements in `results/vlp{17,100}/`, normalises everything into the
benchmarks schema, and writes the two reference CSVs.

Designed to be run once at initial setup. Re-running overwrites the files;
new approved runs / new hand measurements should be appended via dedicated
add scripts (or, for now, by hand).

Usage:
    uv run python -m analysis.seed_benchmarks
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS   = REPO_ROOT / "results"
BENCH_VLP = REPO_ROOT / "benchmarks" / "vlp"

SCRIPT_VERSION = "vlp_measure_v2@2.0(wall=0.75;contrast=0.05;cv=0.2;prom=0.15;smooth=0.7;end=25.0)"


# ── Reference runs (per-image script outputs we trust) ────────────────────

def _per_image_summary(csv_path: Path, batch_id: str, sample_subtype: str,
                       approver: str, notes: str) -> pd.DataFrame:
    """Aggregate a vlp_measurements.csv into one row per image."""
    df = pd.read_csv(csv_path)
    rows = []
    for fname, sub in df.groupby("file", sort=True):
        rel = sub[sub["is_reliable"]]
        n_gold = int(len(sub))
        n_wall = int(sub["capsid_diameter_nm"].notna().sum())
        n_rel  = int(len(rel))
        wcv    = sub["wall_radius_cv"].dropna()
        rows.append({
            "batch_id":               batch_id,
            "filename":               str(fname),
            "script_version":         SCRIPT_VERSION,
            "sample_subtype":         sample_subtype,
            "n_gold":                 n_gold,
            "n_wall_fit":             n_wall,
            "n_reliable":             n_rel,
            "wall_fit_success_rate":  (n_wall / n_gold) if n_gold else np.nan,
            "reliable_rate":          (n_rel  / n_gold) if n_gold else np.nan,
            "drop_rate":              ((n_gold - n_rel) / n_gold) if n_gold else np.nan,
            "gold_mean_nm":           float(rel["gold_diameter_nm"].mean())   if n_rel else np.nan,
            "gold_std_nm":            float(rel["gold_diameter_nm"].std())    if n_rel >= 2 else np.nan,
            "gold_median_nm":         float(rel["gold_diameter_nm"].median()) if n_rel else np.nan,
            "capsid_mean_nm":         float(rel["capsid_diameter_nm"].mean())   if n_rel else np.nan,
            "capsid_std_nm":          float(rel["capsid_diameter_nm"].std())    if n_rel >= 2 else np.nan,
            "capsid_median_nm":       float(rel["capsid_diameter_nm"].median()) if n_rel else np.nan,
            "median_wall_cv":         float(wcv.median())                     if len(wcv)     else np.nan,
            "iqr_wall_cv":            float(wcv.quantile(0.75) - wcv.quantile(0.25)) if len(wcv) >= 2 else np.nan,
            "approved_date":          _dt.date.today().isoformat(),
            "approver":               approver,
            "notes":                  notes,
        })
    return pd.DataFrame(rows)


# ── Reference hand (per-particle hand measurements) ───────────────────────

def _safe_relpath(csv_path: Path) -> str:
    """Return path relative to REPO_ROOT when possible, else absolute."""
    p = csv_path.resolve()
    try:
        return str(p.relative_to(REPO_ROOT))
    except ValueError:
        return str(p)


def _parse_paired_hand(csv_path: Path, length_unit: str, batch_id: str,
                       scientist: str, measure_date: str,
                       notes: str = "") -> pd.DataFrame:
    """
    Hand CSVs from ImageJ where rows alternate gold (odd) / capsid (even),
    so 2N rows = N paired particles. `length_unit` is "um" or "nm".
    """
    raw = pd.read_csv(csv_path)
    if "Length" not in raw.columns:
        raise ValueError(f"{csv_path} has no 'Length' column")
    lengths = raw["Length"].to_numpy()
    if length_unit == "um":
        lengths_nm = lengths * 1000.0
    elif length_unit == "nm":
        lengths_nm = lengths
    else:
        raise ValueError(f"unknown length_unit: {length_unit!r}")

    if len(lengths_nm) % 2 != 0:
        raise ValueError(f"{csv_path} has odd row count ({len(lengths_nm)}); expected paired (gold,capsid)")

    rows = []
    for i in range(0, len(lengths_nm), 2):
        rows.append({
            "batch_id":                   batch_id,
            "image_filename":             None,
            "particle_idx":               i // 2 + 1,
            "hand_gold_diameter_nm":      float(lengths_nm[i]),
            "hand_capsid_diameter_nm":    float(lengths_nm[i + 1]),
            "scientist":                  scientist,
            "measure_date":               measure_date,
            "source_file":                _safe_relpath(csv_path),
            "notes":                      notes,
        })
    return pd.DataFrame(rows)


def _parse_single_hand(csv_path: Path, length_unit: str, batch_id: str,
                       image_filename: str, scientist: str, measure_date: str,
                       notes: str = "") -> pd.DataFrame:
    """Hand CSV where every row is a capsid measurement on the named image."""
    raw = pd.read_csv(csv_path)
    if "Length" not in raw.columns:
        raise ValueError(f"{csv_path} has no 'Length' column")
    lengths = raw["Length"].to_numpy()
    lengths_nm = lengths * 1000.0 if length_unit == "um" else lengths

    rows = []
    for idx, L in enumerate(lengths_nm, start=1):
        rows.append({
            "batch_id":                   batch_id,
            "image_filename":             image_filename,
            "particle_idx":               idx,
            "hand_gold_diameter_nm":      np.nan,  # capsid-only measurement
            "hand_capsid_diameter_nm":    float(L),
            "scientist":                  scientist,
            "measure_date":               measure_date,
            "source_file":                _safe_relpath(csv_path),
            "notes":                      notes,
        })
    return pd.DataFrame(rows)


# ── Main seed routine ──────────────────────────────────────────────────────

def main() -> None:
    BENCH_VLP.mkdir(parents=True, exist_ok=True)

    # 1) Reference runs ────────────────────────────────────────────────────
    runs = pd.concat([
        _per_image_summary(
            csv_path       = RESULTS / "vlp17_v2"  / "vlp_measurements.csv",
            batch_id       = "VLP17_v2_initial",
            sample_subtype = "VLP17",
            approver       = "Lily",
            notes          = "Initial v2 run on VLP17 set; matches hand-measurement aggregate within 0.05 nm. Includes the under-stained 0001-0006 subset (biased ~0.5 nm low) — keep so eval can flag similar quality patterns in future runs.",
        ),
        _per_image_summary(
            csv_path       = RESULTS / "vlp20_v2"  / "vlp_measurements.csv",
            batch_id       = "VLP20_v2_initial",
            sample_subtype = "VLP20",
            approver       = "Lily",
            notes          = "Initial v2 run on VLP20 set. No hand measurements yet for this subtype.",
        ),
        _per_image_summary(
            csv_path       = RESULTS / "vlp100_v2" / "vlp_measurements.csv",
            batch_id       = "VLP_100_v2_initial",
            sample_subtype = "VLP_100",
            approver       = "Lily",
            notes          = "Initial v2 run on VLP_100 set. Reads ~1 nm over hand mean (28.4 hand vs 29.3 script) — possible per-sample tuning target, kept as-is for now.",
        ),
    ], ignore_index=True)

    runs_path = BENCH_VLP / "reference_runs.csv"
    runs.to_csv(runs_path, index=False)
    print(f"reference_runs.csv → {runs_path}  ({len(runs)} rows)")

    # 2) Reference hand ────────────────────────────────────────────────────
    hand = pd.concat([
        _parse_paired_hand(
            csv_path     = RESULTS / "vlp17"  / "VLP_17_hand_measurements.csv",
            length_unit  = "um",
            batch_id     = "VLP17_v2_initial",
            scientist    = "Lily",
            measure_date = "2026-05-01",
            notes        = "Aggregate VLP17 hand measurements (paired gold+capsid, ImageJ). Length in µm in source.",
        ),
        _parse_paired_hand(
            csv_path     = RESULTS / "vlp100" / "VLP_100_hand_measurements.csv",
            length_unit  = "nm",
            batch_id     = "VLP_100_v2_initial",
            scientist    = "Lily",
            measure_date = "2026-05-01",
            notes        = "Aggregate VLP_100 hand measurements (paired gold+capsid, ImageJ). Length in nm in source.",
        ),
        # Per-image diagnostics (vlp17_003_hand.csv, vlp17_0010_hand.csv) intentionally
        # excluded from the reference set — those are debugging measurements, not
        # broad sample-batch ground truth. They live in results/vlp17/ for
        # reproducibility but shouldn't drive eval thresholds.
    ], ignore_index=True)

    hand_path = BENCH_VLP / "reference_hand.csv"
    hand.to_csv(hand_path, index=False)
    print(f"reference_hand.csv → {hand_path}  ({len(hand)} rows)")

    # 3) Quick summary so the seed can be sanity-checked ───────────────────
    print("\n── seeded reference_runs.csv ──")
    print(runs.groupby("batch_id").agg(
        n_images=("filename", "count"),
        capsid_mean_nm=("capsid_mean_nm", "mean"),
        reliable_rate=("reliable_rate", "mean"),
    ).round(3).to_string())

    print("\n── seeded reference_hand.csv ──")
    summary_rows = []
    for bid, sub in hand.groupby("batch_id"):
        per_image = sub.groupby(sub["image_filename"].fillna("(aggregate)"))
        for img, ssub in per_image:
            summary_rows.append({
                "batch_id":     bid,
                "image":        img,
                "n_particles":  len(ssub),
                "gold_mean_nm": float(ssub["hand_gold_diameter_nm"].mean()) if ssub["hand_gold_diameter_nm"].notna().any() else np.nan,
                "capsid_mean_nm": float(ssub["hand_capsid_diameter_nm"].mean()),
            })
    print(pd.DataFrame(summary_rows).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
