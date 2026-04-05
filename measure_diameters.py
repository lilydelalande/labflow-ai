"""
measure_diameters.py — TEM particle size analysis

Reads .dm3 / .dm4 Gatan Digital Micrograph files, detects particles using
Difference-of-Gaussians blob detection, measures equivalent circular diameters
in nm, and outputs:
  - overlay image showing detections (PNG)
  - size histogram (PNG)
  - CSV with per-particle measurements

Usage:
    uv run python measure_diameters.py <image_dir> [options]

Examples:
    uv run python measure_diameters.py "images/VLPs for machine learning project"
    uv run python measure_diameters.py "images/" --pattern "VLP17_*" --out results/
    uv run python measure_diameters.py "images/" --min-nm 20 --max-nm 200
    uv run python measure_diameters.py "images/" --threshold 0.02
"""

import argparse
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import ncempy.io as nio
import numpy as np
import pandas as pd
from skimage import exposure
from skimage.feature import blob_dog


# ── Image loading ─────────────────────────────────────────────────────────────

def load_dm(path: Path) -> tuple[np.ndarray, float]:
    """Return (image as float32, nm_per_pixel)."""
    d = nio.read(str(path))
    img = d["data"].astype(np.float32)
    nm_per_px = float(d["pixelSize"][0])
    return img, nm_per_px


def normalise(img: np.ndarray) -> np.ndarray:
    """Percentile-stretch to [0,1] then CLAHE for local contrast enhancement."""
    lo, hi = np.percentile(img, [0.5, 99.5])
    img = np.clip((img - lo) / (hi - lo), 0, 1)
    return exposure.equalize_adapthist(img, clip_limit=0.02)


# ── Detection ─────────────────────────────────────────────────────────────────

def detect_blobs(
    img_norm: np.ndarray,
    nm_per_px: float,
    threshold: float,
    min_nm: float | None,
    max_nm: float | None,
    dark_on_bright: bool = True,
) -> pd.DataFrame:
    """
    Detect circular particles using Difference-of-Gaussians blob detection.

    DoG searches across multiple scales simultaneously — no need to know particle
    size in advance. Each detected blob gives a centroid and an approximate radius.

    Parameters
    ----------
    img_norm      : float image in [0, 1]
    nm_per_px     : physical pixel size
    threshold     : detection sensitivity (lower = more detections, more noise).
                    Start at 0.05 and decrease toward 0.01 if particles are missed.
    min_nm/max_nm : optional size guard rails (nm). Leave None to see everything.
    dark_on_bright: True for negative stain (dark particles on bright background),
                    False for cryo / bright-on-dark images.
    """
    # DoG finds bright blobs; invert if particles are dark
    img_for_blob = (1.0 - img_norm) if dark_on_bright else img_norm

    # Sigma (px) relates to blob radius: radius = sigma * sqrt(2)
    # Search from 2 nm radius up to 500 nm radius by default
    min_sigma_px = max(1.0, (min_nm / 2 if min_nm else 2.0) / nm_per_px / np.sqrt(2))
    max_sigma_px = ((max_nm / 2 if max_nm else 500.0) / nm_per_px / np.sqrt(2))

    blobs = blob_dog(
        img_for_blob,
        min_sigma=min_sigma_px,
        max_sigma=max_sigma_px,
        sigma_ratio=1.6,
        threshold=threshold,
    )

    if len(blobs) == 0:
        return pd.DataFrame(columns=["centroid-0", "centroid-1", "diameter_nm"])

    df = pd.DataFrame(blobs, columns=["centroid-0", "centroid-1", "sigma"])
    df["diameter_nm"] = 2 * df["sigma"] * np.sqrt(2) * nm_per_px

    if min_nm is not None:
        df = df[df["diameter_nm"] >= min_nm]
    if max_nm is not None:
        df = df[df["diameter_nm"] <= max_nm]

    return df.reset_index(drop=True)


# ── Output ────────────────────────────────────────────────────────────────────

def save_overlay(
    img_norm: np.ndarray,
    df: pd.DataFrame,
    nm_per_px: float,
    out_path: Path,
    title: str = "",
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    axes[0].imshow(img_norm, cmap="gray", interpolation="none")
    axes[0].set_title("Image")
    axes[0].axis("off")

    axes[1].imshow(img_norm, cmap="gray", interpolation="none")
    for _, row in df.iterrows():
        cy, cx = row["centroid-0"], row["centroid-1"]
        r_px = row["diameter_nm"] / 2 / nm_per_px
        circle = mpatches.Circle(
            (cx, cy), r_px, linewidth=0.8, edgecolor="cyan", facecolor="none"
        )
        axes[1].add_patch(circle)
    axes[1].set_title(f"Detections (n={len(df)})")
    axes[1].axis("off")

    if title:
        fig.suptitle(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_histogram(df_all: pd.DataFrame, out_path: Path) -> None:
    diameters = df_all["diameter_nm"]
    mean_d = diameters.mean()
    std_d = diameters.std()

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(diameters, bins=40, edgecolor="white", color="steelblue")
    ax.axvline(mean_d, color="tomato", linewidth=2,
               label=f"mean = {mean_d:.1f} ± {std_d:.1f} nm  (n={len(diameters)})")
    ax.set_xlabel("Equivalent diameter (nm)")
    ax.set_ylabel("Count")
    ax.set_title("Particle size distribution")
    ax.legend()
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Histogram  → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("image_dir", type=Path,
                        help="Folder containing .dm3 / .dm4 files")
    parser.add_argument("--pattern", default="*",
                        help="Filename glob, e.g. 'VLP17_*' (default: all dm files)")
    parser.add_argument("--out", type=Path, default=Path("results"),
                        help="Output directory (default: results/)")
    parser.add_argument("--min-nm", type=float, default=None,
                        help="Optional: exclude particles smaller than this (nm)")
    parser.add_argument("--max-nm", type=float, default=None,
                        help="Optional: exclude particles larger than this (nm)")
    parser.add_argument("--threshold", type=float, default=0.05,
                        help="Detection threshold (default 0.05; lower = more sensitive)")
    parser.add_argument("--bright-on-dark", action="store_true",
                        help="Particles are bright on dark background (e.g. cryo-TEM)")
    args = parser.parse_args()

    files = sorted(
        p
        for ext in (".dm3", ".dm4")
        for p in args.image_dir.glob(f"{args.pattern}{ext}")
    )
    if not files:
        print(f"No .dm3 / .dm4 files found in {args.image_dir}", file=sys.stderr)
        sys.exit(1)

    args.out.mkdir(parents=True, exist_ok=True)
    overlay_dir = args.out / "overlays"
    overlay_dir.mkdir(exist_ok=True)

    all_rows = []

    for path in files:
        print(f"\n{path.name}")
        img, nm_per_px = load_dm(path)
        print(f"  {img.shape[0]}×{img.shape[1]} px  |  {nm_per_px:.4f} nm/px")

        img_norm = normalise(img)

        df = detect_blobs(
            img_norm, nm_per_px,
            threshold=args.threshold,
            min_nm=args.min_nm,
            max_nm=args.max_nm,
            dark_on_bright=not args.bright_on_dark,
        )
        df.insert(0, "file", path.name)
        all_rows.append(df)

        n = len(df)
        if n:
            mean_d = df["diameter_nm"].mean()
            std_d = df["diameter_nm"].std()
            print(f"  {n} particles  |  mean {mean_d:.1f} ± {std_d:.1f} nm")
        else:
            print("  No particles detected — try lowering --threshold")

        overlay_path = overlay_dir / (path.stem + "_overlay.png")
        save_overlay(img_norm, df, nm_per_px, overlay_path, title=path.name)
        print(f"  Overlay    → {overlay_path}")

    if not all_rows or all(len(r) == 0 for r in all_rows):
        print("\nNo particles detected across any file.")
        return

    results = pd.concat(all_rows, ignore_index=True)

    csv_path = args.out / "diameters.csv"
    results.to_csv(csv_path, index=False)
    print(f"\n  CSV        → {csv_path}")

    hist_path = args.out / "size_histogram.png"
    save_histogram(results, hist_path)

    print("\n── Summary ──────────────────────────────────────────────")
    print(results["diameter_nm"].describe().rename("diameter (nm)").to_string())

    per_file = (
        results.groupby("file")["diameter_nm"]
        .agg(n="count", mean="mean", std="std", median="median")
        .round(1)
    )
    print("\nPer-file:")
    print(per_file.to_string())


if __name__ == "__main__":
    main()
