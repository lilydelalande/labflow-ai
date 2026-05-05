"""
plot_capsid_groups.py — Per-image gold + capsid strip plot

Reads a vlp_measurements.csv and produces a two-panel strip plot showing
per-image median gold NP and capsid diameters side by side, to distinguish
imaging artifacts (both shift) from real biology (capsid shifts, gold stable).

Usage:
    uv run python plot_capsid_groups.py results/vlp17/vlp_measurements.csv
    uv run python plot_capsid_groups.py results/vlp17/vlp_measurements.csv --out results/vlp17/
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("csv", type=Path)
    parser.add_argument("--out", type=Path, default=None,
                        help="Output directory (default: same folder as CSV)")
    args = parser.parse_args()

    if not args.csv.exists():
        print(f"Not found: {args.csv}", file=sys.stderr)
        sys.exit(1)

    out_dir = args.out or args.csv.parent
    df = pd.read_csv(args.csv)

    if "is_reliable" in df.columns:
        df = df[df["is_reliable"]]
    df = df.dropna(subset=["capsid_diameter_nm"])

    files = sorted(df["file"].unique())

    per_image = (
        df.groupby("file")
        .agg(
            capsid_median=("capsid_diameter_nm", "median"),
            capsid_std=("capsid_diameter_nm", "std"),
            gold_median=("gold_diameter_nm", "median"),
            gold_std=("gold_diameter_nm", "std"),
            n=("capsid_diameter_nm", "count"),
        )
        .reindex(files)
        .reset_index()
    )
    per_image["file_i"] = range(len(per_image))
    xlabels = [Path(f).stem for f in per_image["file"]]

    # ── Strip plot ────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    for ax, col, std_col, color, ecolor, ylabel, title in [
        (axes[0], "gold_median", "gold_std", "goldenrod", "wheat",
         "Gold NP diameter (nm)", "Per-image median gold NP diameter ± std"),
        (axes[1], "capsid_median", "capsid_std", "steelblue", "lightsteelblue",
         "Capsid diameter (nm)", "Per-image median capsid diameter ± std"),
    ]:
        ax.errorbar(
            per_image["file_i"], per_image[col],
            yerr=per_image[std_col],
            fmt="o", markersize=6, color=color, ecolor=ecolor,
            capsize=3, linewidth=1,
        )
        ax.axhline(per_image[col].median(), color="tomato", linestyle="--",
                   linewidth=1, label=f"overall median = {per_image[col].median():.1f} nm")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)

        # Y-axis: zoom to the main cluster (5th–95th percentile ± 20%)
        lo = per_image[col].quantile(0.05)
        hi = per_image[col].quantile(0.95)
        pad = (hi - lo) * 0.6 + 0.5
        ax.set_ylim(lo - pad, hi + pad)

    axes[1].set_xticks(per_image["file_i"])
    axes[1].set_xticklabels(xlabels, rotation=45, ha="right", fontsize=7)

    plt.tight_layout()
    out_path = out_dir / "capsid_groups.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Strip plot → {out_path}")

    print("\n── Per-image medians ────────────────────────────────────")
    print(per_image[["file", "gold_median", "capsid_median", "capsid_std", "n"]].to_string(index=False))


if __name__ == "__main__":
    main()
