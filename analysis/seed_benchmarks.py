"""
seed_benchmarks.py — Initial seed of benchmarks/<sample_type>/{reference_runs.csv,reference_hand.csv}.

Reads the existing v2 outputs in `results/vlp{17,20,100}_v2/` and the hand
measurements in `results/vlp{17,100}/`, normalises everything into the
benchmarks schema, and writes the two reference CSVs.

Schema (v2 — drops the old batch_id and sample_subtype columns):

  reference_runs.csv:
    sample_name, filename, script_version,
    n_gold, n_wall_fit, n_reliable,
    wall_fit_success_rate, reliable_rate, drop_rate,
    gold_{mean,std,median}_nm, capsid_{mean,std,median}_nm,
    median_wall_cv, iqr_wall_cv,
    approved_date, approver, notes

  reference_hand.csv:
    sample_name, particle_idx,
    hand_gold_diameter_nm, hand_capsid_diameter_nm,
    scientist, measure_date, source_file, notes

`sample_name` is the run's input folder name (the same string `run()` writes
into its result dict). No "subtype" column anywhere — eval pools across the
whole sample_type, judging quality (size-invariant rates and CVs) rather
than absolute size.

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
BENCH_BMV = REPO_ROOT / "benchmarks" / "bmv"

SCRIPT_VERSION_VLP = "vlp_measure_v2@2.0(wall=0.75;contrast=0.05;cv=0.2;prom=0.15;smooth=0.7;end=25.0)"
SCRIPT_VERSION_BMV = "bmv_measure@1.0"


# ── Reference runs (per-image script outputs we trust) ────────────────────

def _per_image_summary(csv_path: Path, sample_name: str,
                       approver: str, notes: str,
                       script_version: str = SCRIPT_VERSION_VLP) -> pd.DataFrame:
    """Aggregate a measurements CSV into one row per image. Works for both VLP
    (gold + capsid) and BMV (capsid only — gold columns will be NaN)."""
    df = pd.read_csv(csv_path)
    has_gold = "gold_diameter_nm" in df.columns
    rows = []
    for fname, sub in df.groupby("file", sort=True):
        rel = sub[sub["is_reliable"]]
        n_gold = int(len(sub))
        n_wall = int(sub["capsid_diameter_nm"].notna().sum())
        n_rel  = int(len(rel))
        wcv    = sub["wall_radius_cv"].dropna()
        rows.append({
            "sample_name":            sample_name,
            "filename":               str(fname),
            "script_version":         script_version,
            "n_gold":                 n_gold,
            "n_wall_fit":             n_wall,
            "n_reliable":             n_rel,
            "wall_fit_success_rate":  (n_wall / n_gold) if n_gold else np.nan,
            "reliable_rate":          (n_rel  / n_gold) if n_gold else np.nan,
            "drop_rate":              ((n_gold - n_rel) / n_gold) if n_gold else np.nan,
            "gold_mean_nm":           float(rel["gold_diameter_nm"].mean())   if (has_gold and n_rel) else np.nan,
            "gold_std_nm":            float(rel["gold_diameter_nm"].std())    if (has_gold and n_rel >= 2) else np.nan,
            "gold_median_nm":         float(rel["gold_diameter_nm"].median()) if (has_gold and n_rel) else np.nan,
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


def _parse_paired_hand(csv_path: Path, length_unit: str, sample_name: str,
                       scientist: str, measure_date: str,
                       notes: str = "") -> pd.DataFrame:
    """ImageJ CSV with rows alternating gold (odd) / capsid (even). 2N rows = N particles."""
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
            "sample_name":             sample_name,
            "particle_idx":            i // 2 + 1,
            "hand_gold_diameter_nm":   float(lengths_nm[i]),
            "hand_capsid_diameter_nm": float(lengths_nm[i + 1]),
            "scientist":               scientist,
            "measure_date":            measure_date,
            "source_file":             _safe_relpath(csv_path),
            "notes":                   notes,
        })
    return pd.DataFrame(rows)


def _parse_single_hand(csv_path: Path, sample_name: str,
                       scientist: str, measure_date: str,
                       notes: str = "") -> pd.DataFrame:
    """Single-column ImageJ CSV (one capsid measurement per row, no gold).
    Per-row unit detection — values < 1.0 are µm, otherwise nm. Handles
    mixed-unit exports where the image scale changed mid-measurement."""
    raw = pd.read_csv(csv_path)
    if "Length" not in raw.columns:
        raise ValueError(f"{csv_path} has no 'Length' column")
    lengths = raw["Length"].to_numpy(dtype=float)
    lengths_nm = np.where(lengths < 1.0, lengths * 1000.0, lengths)

    rows = []
    for i, L in enumerate(lengths_nm, start=1):
        rows.append({
            "sample_name":             sample_name,
            "particle_idx":            i,
            "hand_gold_diameter_nm":   np.nan,   # no gold in single-column hand
            "hand_capsid_diameter_nm": float(L),
            "scientist":               scientist,
            "measure_date":            measure_date,
            "source_file":             _safe_relpath(csv_path),
            "notes":                   notes,
        })
    return pd.DataFrame(rows)


# ── Main seed routine ──────────────────────────────────────────────────────

def main() -> None:
    BENCH_VLP.mkdir(parents=True, exist_ok=True)
    BENCH_BMV.mkdir(parents=True, exist_ok=True)

    # ── VLP ─────────────────────────────────────────────────────────────────
    vlp_runs = pd.concat([
        _per_image_summary(
            csv_path    = RESULTS / "vlp17_v2"  / "vlp_measurements.csv",
            sample_name = "VLP17_v2_initial",
            approver    = "Lily",
            notes       = "Initial v2 run on VLP17 set; matches hand-measurement aggregate within 0.05 nm. Includes the under-stained 0001-0006 subset (biased ~0.5 nm low) — keep so eval can flag similar quality patterns in future runs.",
        ),
        _per_image_summary(
            csv_path    = RESULTS / "vlp20_v2"  / "vlp_measurements.csv",
            sample_name = "VLP20_v2_initial",
            approver    = "Lily",
            notes       = "Initial v2 run on VLP20 set. No hand measurements yet.",
        ),
        _per_image_summary(
            csv_path    = RESULTS / "vlp100_v2" / "vlp_measurements.csv",
            sample_name = "VLP_100_v2_initial",
            approver    = "Lily",
            notes       = "Initial v2 run on VLP_100 set. Reads ~1 nm over hand mean (28.4 hand vs 29.3 script) — possible per-sample tuning target, kept as-is for now.",
        ),
    ], ignore_index=True)

    vlp_runs_path = BENCH_VLP / "reference_runs.csv"
    vlp_runs.to_csv(vlp_runs_path, index=False)
    print(f"VLP reference_runs.csv → {vlp_runs_path}  ({len(vlp_runs)} rows)")

    vlp_hand = pd.concat([
        _parse_paired_hand(
            csv_path     = RESULTS / "vlp17"  / "VLP_17_hand_measurements.csv",
            length_unit  = "um",
            sample_name  = "VLP17_v2_initial",
            scientist    = "Lily",
            measure_date = "2026-05-01",
            notes        = "Aggregate VLP17 hand measurements (paired gold+capsid, ImageJ). Length in µm in source.",
        ),
        _parse_paired_hand(
            csv_path     = RESULTS / "vlp100" / "VLP_100_hand_measurements.csv",
            length_unit  = "nm",
            sample_name  = "VLP_100_v2_initial",
            scientist    = "Lily",
            measure_date = "2026-05-01",
            notes        = "Aggregate VLP_100 hand measurements (paired gold+capsid, ImageJ). Length in nm in source.",
        ),
    ], ignore_index=True)

    vlp_hand_path = BENCH_VLP / "reference_hand.csv"
    vlp_hand.to_csv(vlp_hand_path, index=False)
    print(f"VLP reference_hand.csv → {vlp_hand_path}  ({len(vlp_hand)} rows)")

    # ── BMV ─────────────────────────────────────────────────────────────────
    bmv_runs = _per_image_summary(
        csv_path       = RESULTS / "bmv" / "bmv_measurements.csv",
        sample_name    = "BMV_initial",
        approver       = "Lily",
        notes          = "Initial bmv_measure run; capsid mean ~28.6 nm matches literature.",
        script_version = SCRIPT_VERSION_BMV,
    )
    bmv_runs_path = BENCH_BMV / "reference_runs.csv"
    bmv_runs.to_csv(bmv_runs_path, index=False)
    print(f"BMV reference_runs.csv → {bmv_runs_path}  ({len(bmv_runs)} rows)")

    bmv_hand_csv = RESULTS / "bmv" / "hand" / "bmv_hand_measurements.csv"
    if bmv_hand_csv.exists():
        bmv_hand = _parse_single_hand(
            csv_path     = bmv_hand_csv,
            sample_name  = "BMV_initial",
            scientist    = "Lily",
            measure_date = "2026-05-06",
            notes        = "Aggregate BMV hand measurements (capsid only — no gold). Mixed nm/µm units in source; per-row unit detection handles it.",
        )
        bmv_hand_path = BENCH_BMV / "reference_hand.csv"
        bmv_hand.to_csv(bmv_hand_path, index=False)
        print(f"BMV reference_hand.csv → {bmv_hand_path}  ({len(bmv_hand)} rows)")

    print("\n── seeded VLP reference_runs.csv ──")
    print(vlp_runs.groupby("sample_name").agg(
        n_images=("filename", "count"),
        capsid_mean_nm=("capsid_mean_nm", "mean"),
        reliable_rate=("reliable_rate", "mean"),
    ).round(3).to_string())

    print("\n── seeded VLP reference_hand.csv ──")
    print(vlp_hand.groupby("sample_name").agg(
        n_particles=("particle_idx", "count"),
        gold_mean_nm=("hand_gold_diameter_nm", "mean"),
        capsid_mean_nm=("hand_capsid_diameter_nm", "mean"),
    ).round(2).to_string())

    print("\n── seeded BMV reference_runs.csv ──")
    print(bmv_runs.groupby("sample_name").agg(
        n_images=("filename", "count"),
        capsid_mean_nm=("capsid_mean_nm", "mean"),
        reliable_rate=("reliable_rate", "mean"),
    ).round(3).to_string())
    if bmv_hand_csv.exists():
        print("\n── seeded BMV reference_hand.csv ──")
        print(bmv_hand.groupby("sample_name").agg(
            n_particles=("particle_idx", "count"),
            capsid_mean_nm=("hand_capsid_diameter_nm", "mean"),
            capsid_std_nm=("hand_capsid_diameter_nm", "std"),
        ).round(2).to_string())


if __name__ == "__main__":
    main()
