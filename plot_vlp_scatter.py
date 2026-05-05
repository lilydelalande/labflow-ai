"""
plot_vlp_scatter.py — Combined scatter + 2D histogram across multiple VLP samples

Reads one or more vlp_measurements.csv files, infers sample group from the
filename column, and produces:
  - combined_scatter.png   : scatter plot, one colour per sample
  - combined_hist2d.png    : 2D histogram (heatmap) per sample in subplots

Usage:
    uv run python plot_vlp_scatter.py results/vlp_measurements.csv
    uv run python plot_vlp_scatter.py vlp17/vlp_measurements.csv vlp20/vlp_measurements.csv vlp100/vlp_measurements.csv
    uv run python plot_vlp_scatter.py results/vlp_measurements.csv --out-dir results/
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


COLORS = ["steelblue", "tomato", "seagreen", "mediumpurple", "darkorange", "hotpink"]


def infer_sample(filename: str) -> str:
    """
    Extract sample name from a TEM image filename.
    E.g. VLP17_0001.dm4 → VLP17, VLP_100_06_14_0001.dm3 → VLP_100
    Matches the leading letters/underscores + digits block.
    """
    stem = Path(filename).stem
    m = re.match(r"^([A-Za-z_]+\d*)", stem)
    return m.group(1) if m else stem


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("csvs", nargs="+", type=Path,
                        help="One or more vlp_measurements.csv files")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory (default: same folder as first CSV)")
    parser.add_argument("--reliable-only", action="store_true", default=True,
                        help="Plot only reliable detections (default: True)")
    parser.add_argument("--combined-hist2d", action="store_true",
                        help="Plot all samples on one 2D histogram instead of per-sample subplots")
    args = parser.parse_args()

    frames = []
    for p in args.csvs:
        if not p.exists():
            print(f"Not found: {p}", file=sys.stderr)
            sys.exit(1)
        df = pd.read_csv(p)
        frames.append(df)

    results = pd.concat(frames, ignore_index=True)
    results["sample"] = results["file"].apply(infer_sample)

    # Filter to reliable only
    if "is_reliable" in results.columns:
        plot_data = results[results["is_reliable"]].copy()
        n_flagged = (~results["is_reliable"]).sum()
    else:
        plot_data = results.copy()
        n_flagged = 0

    plot_data = plot_data.dropna(subset=["capsid_diameter_nm", "gold_diameter_nm"])

    samples = sorted(plot_data["sample"].unique())

    fig, ax = plt.subplots(figsize=(8, 7))

    for i, sample in enumerate(samples):
        color = COLORS[i % len(COLORS)]
        sub = plot_data[plot_data["sample"] == sample]
        ax.scatter(
            sub["gold_diameter_nm"],
            sub["capsid_diameter_nm"],
            alpha=0.45,
            s=18,
            color=color,
            edgecolors="none",
            label=f"{sample}  (n={len(sub)}, gold {sub['gold_diameter_nm'].mean():.1f}±{sub['gold_diameter_nm'].std():.1f} nm, "
                  f"capsid {sub['capsid_diameter_nm'].mean():.1f}±{sub['capsid_diameter_nm'].std():.1f} nm)",
        )
        # Mean crosshair per sample
        ax.axvline(sub["gold_diameter_nm"].mean(), color=color, linewidth=0.8,
                   linestyle="--", alpha=0.6)
        ax.axhline(sub["capsid_diameter_nm"].mean(), color=color, linewidth=0.8,
                   linestyle="--", alpha=0.6)

    n_total = len(plot_data)
    flag_note = f"  ({n_flagged} flagged excluded)" if n_flagged else ""
    title = f"VLP capsid vs gold NP diameter — {n_total} particles{flag_note}"
    ax.set_xlabel("Gold NP diameter (nm)")
    ax.set_ylabel("VLP capsid diameter (nm)")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()

    out_dir = args.out_dir or args.csvs[0].parent
    scatter_path = out_dir / "combined_scatter.png"
    fig.savefig(scatter_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Scatter   → {scatter_path}")

    # ── 2D histogram ─────────────────────────────────────────────────────────
    x_all = plot_data["gold_diameter_nm"]
    y_all = plot_data["capsid_diameter_nm"]
    x_lim = (np.floor(x_all.min()) - 1, np.ceil(x_all.max()) + 1)
    y_lim = (np.floor(y_all.min()) - 1, np.ceil(y_all.max()) + 1)

    if args.combined_hist2d:
        # Normalize each sample to its own density so high-n samples don't
        # swamp low-n ones. Sum the normalized histograms before plotting.
        bins_x = np.linspace(x_lim[0], x_lim[1], 31)
        bins_y = np.linspace(y_lim[0], y_lim[1], 31)
        combined = np.zeros((30, 30))
        for sample in samples:
            sub = plot_data[plot_data["sample"] == sample]
            h, _, _ = np.histogram2d(
                sub["gold_diameter_nm"], sub["capsid_diameter_nm"],
                bins=[bins_x, bins_y],
            )
            if h.max() > 0:
                combined += h / h.max()  # normalise to peak = 1 per sample

        fig2, ax2 = plt.subplots(figsize=(7, 6))
        im = ax2.pcolormesh(bins_x, bins_y, combined.T, cmap="viridis")
        fig2.colorbar(im, ax=ax2, label="normalised intensity (peak=1 per sample)")
        ax2.set_xlabel("Gold NP diameter (nm)")
        ax2.set_ylabel("Capsid diameter (nm)")
        ax2.set_title(f"All samples — normalised 2D histogram{flag_note}\n" +
                      "  ".join(f"{s} n={len(plot_data[plot_data['sample']==s])}" for s in samples))
        plt.tight_layout()
        hist2d_path = out_dir / "combined_hist2d.png"
        fig2.savefig(hist2d_path, dpi=150, bbox_inches="tight")
        plt.close(fig2)
        print(f"2D hist   → {hist2d_path}")
    else:
        cols = min(len(samples), 3)
        rows = (len(samples) + cols - 1) // cols
        fig2, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows), squeeze=False)

        for i, sample in enumerate(samples):
            ax2 = axes[i // cols][i % cols]
            sub = plot_data[plot_data["sample"] == sample]
            h = ax2.hist2d(
                sub["gold_diameter_nm"],
                sub["capsid_diameter_nm"],
                bins=30,
                range=[x_lim, y_lim],
                cmap="viridis",
            )
            fig2.colorbar(h[3], ax=ax2, label="count")
            ax2.set_xlabel("Gold NP diameter (nm)")
            ax2.set_ylabel("Capsid diameter (nm)")
            ax2.set_title(
                f"{sample}  (n={len(sub)})\n"
                f"gold {sub['gold_diameter_nm'].mean():.1f}±{sub['gold_diameter_nm'].std():.1f} nm  "
                f"capsid {sub['capsid_diameter_nm'].mean():.1f}±{sub['capsid_diameter_nm'].std():.1f} nm"
            )

        for j in range(i + 1, rows * cols):
            axes[j // cols][j % cols].set_visible(False)

        fig2.suptitle(f"VLP capsid vs gold NP — 2D histogram{flag_note}", fontsize=12)
        plt.tight_layout()
        hist2d_path = out_dir / "combined_hist2d.png"
        fig2.savefig(hist2d_path, dpi=150, bbox_inches="tight")
        plt.close(fig2)
        print(f"2D hist   → {hist2d_path}")

    # ── KDE contour plot ─────────────────────────────────────────────────────
    fig3, ax3 = plt.subplots(figsize=(8, 7))

    xx, yy = np.mgrid[x_lim[0]:x_lim[1]:200j, y_lim[0]:y_lim[1]:200j]
    grid_points = np.vstack([xx.ravel(), yy.ravel()])

    for i, sample in enumerate(samples):
        color = COLORS[i % len(COLORS)]
        sub = plot_data[plot_data["sample"] == sample]
        if len(sub) < 4:
            continue
        kde = gaussian_kde(
            np.vstack([sub["gold_diameter_nm"], sub["capsid_diameter_nm"]]),
            bw_method="scott",
        )
        zz = kde(grid_points).reshape(xx.shape)
        zz /= zz.max()  # normalise each sample to peak=1
        ax3.contour(xx, yy, zz, levels=[0.1, 0.3, 0.5, 0.7, 0.9],
                    colors=[color], linewidths=1.2, alpha=0.85)
        ax3.contourf(xx, yy, zz, levels=[0.5, 1.0],
                     colors=[color], alpha=0.15)
        # Mark the peak
        peak = sub[["gold_diameter_nm", "capsid_diameter_nm"]].mean()
        ax3.plot(peak["gold_diameter_nm"], peak["capsid_diameter_nm"],
                 "o", color=color, markersize=6,
                 label=f"{sample} (n={len(sub)})")

    ax3.set_xlabel("Gold NP diameter (nm)")
    ax3.set_ylabel("Capsid diameter (nm)")
    ax3.set_title(f"VLP capsid vs gold NP — KDE contours{flag_note}\n"
                  "(contours at 10/30/50/70/90% of each sample's peak density)")
    ax3.legend(fontsize=9)
    plt.tight_layout()
    kde_path = out_dir / "combined_kde.png"
    fig3.savefig(kde_path, dpi=150, bbox_inches="tight")
    plt.close(fig3)
    print(f"KDE       → {kde_path}")

    # ── Pooled 1D capsid histogram ────────────────────────────────────────────
    fig4, ax4 = plt.subplots(figsize=(9, 4))
    bins = np.arange(
        np.floor(plot_data["capsid_diameter_nm"].min()),
        np.ceil(plot_data["capsid_diameter_nm"].max()) + 1,
        1,
    )
    for i, sample in enumerate(samples):
        sub = plot_data[plot_data["sample"] == sample]
        ax4.hist(sub["capsid_diameter_nm"], bins=bins, alpha=0.55,
                 color=COLORS[i % len(COLORS)], edgecolor="white",
                 linewidth=0.3, label=f"{sample} (n={len(sub)})")
    ax4.set_xlabel("Capsid diameter (nm)")
    ax4.set_ylabel("Count")
    ax4.set_title(f"Pooled capsid diameter distribution — all samples{flag_note}\n"
                  "1 nm bins  |  if T-number polymorphism exists, expect peaks at fixed absolute sizes")
    ax4.legend(fontsize=9)
    plt.tight_layout()
    hist1d_path = out_dir / "combined_capsid_hist.png"
    fig4.savefig(hist1d_path, dpi=150, bbox_inches="tight")
    plt.close(fig4)
    print(f"Capsid hist → {hist1d_path}")

    print("\n── Per-sample summary ───────────────────────────────────")
    for sample in samples:
        sub = plot_data[plot_data["sample"] == sample]
        print(f"\n{sample}  (n={len(sub)})")
        print(f"  gold   {sub['gold_diameter_nm'].mean():.1f} ± {sub['gold_diameter_nm'].std():.1f} nm")
        print(f"  capsid {sub['capsid_diameter_nm'].mean():.1f} ± {sub['capsid_diameter_nm'].std():.1f} nm")


if __name__ == "__main__":
    main()
