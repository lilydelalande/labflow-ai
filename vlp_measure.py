"""
vlp_measure.py — Gold nanoparticle detection in VLP TEM images

Detects gold nanoparticles (dense, near-black circular blobs) in Gatan
.dm3 / .dm4 TEM images and reports their diameters in nm.

Outputs:
  - overlay image with detected NPs circled (PNG)
  - size histogram (PNG)
  - CSV with per-particle measurements

Usage:
    uv run python vlp_measure.py <image_path> [options]

Examples:
    uv run python vlp_measure.py "images/VLPs for machine learning project/VLP17_0001.dm4"
    uv run python vlp_measure.py "images/VLPs for machine learning project" --pattern "VLP17_*"
    uv run python vlp_measure.py "images/" --gold-threshold 0.04 --min-gold-nm 12
"""

import argparse
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import ncempy.io as nio
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import exposure, filters, morphology, measure


# ── Image loading ─────────────────────────────────────────────────────────────

def load_dm(path: Path) -> tuple[np.ndarray, float]:
    """Return (image as float32, nm_per_pixel)."""
    d = nio.read(str(path))
    img = d["data"].astype(np.float32)
    nm_per_px = float(d["pixelSize"][0])
    return img, nm_per_px


def normalise(img: np.ndarray) -> np.ndarray:
    """Percentile-stretch to [0,1] then CLAHE for local contrast."""
    lo, hi = np.percentile(img, [0.5, 99.5])
    img = np.clip((img - lo) / (hi - lo), 0, 1)
    return exposure.equalize_adapthist(img, clip_limit=0.02)


# ── Gold NP detection ─────────────────────────────────────────────────────────

def auto_threshold(img_norm: np.ndarray, gold_threshold: float | None) -> float:
    """
    Find the intensity cutoff for gold NPs automatically per image.

    Uses the valley between the dark gold-NP peak and the brighter background
    in the intensity histogram (multi-Otsu on the dark half of the image).
    Falls back to the user-supplied value if automatic detection fails.
    """
    if gold_threshold is not None:
        return gold_threshold

    # Only look at the dark half of the histogram where gold NPs live
    dark_pixels = img_norm[img_norm < 0.5]
    if len(dark_pixels) < 100:
        return 0.05

    try:
        # Two-class Otsu on the dark pixels finds the valley between noise and gold
        thresh = filters.threshold_otsu(dark_pixels)
        return float(thresh)
    except Exception:
        return 0.05


def detect_gold(
    img_norm: np.ndarray,
    nm_per_px: float,
    threshold: float | None,
    min_gold_nm: float,
    max_gold_nm: float,
) -> tuple[pd.DataFrame, float]:
    """
    Detect gold NPs by intensity thresholding.

    Gold NPs are the darkest feature in the image. The threshold is found
    automatically per image (recommended) or can be set manually.
    Results are filtered to [min_gold_nm, max_gold_nm] to reject noise
    and aggregates.

    Returns (DataFrame: centroid_y, centroid_x, diameter_nm), threshold_used.
    """
    t = auto_threshold(img_norm, threshold)
    binary = img_norm < t
    binary = morphology.remove_small_objects(binary, max_size=3)

    labels = ndi.label(binary)[0]
    props = measure.regionprops_table(
        labels,
        properties=("label", "area", "centroid", "eccentricity", "solidity"),
    )
    df = pd.DataFrame(props)
    if df.empty:
        return pd.DataFrame(columns=["centroid_y", "centroid_x", "diameter_nm"])

    df["diameter_nm"] = 2 * np.sqrt(df["area"] / np.pi) * nm_per_px

    df = df[
        (df["solidity"] > 0.75)
        & (df["eccentricity"] < 0.75)
        & (df["diameter_nm"] >= min_gold_nm)
        & (df["diameter_nm"] <= max_gold_nm)
    ]

    df = df.rename(columns={"centroid-0": "centroid_y", "centroid-1": "centroid_x"})
    return df[["centroid_y", "centroid_x", "diameter_nm"]].reset_index(drop=True), t


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
        cx, cy = row["centroid_x"], row["centroid_y"]
        r_px = row["diameter_nm"] / 2 / nm_per_px
        axes[1].add_patch(mpatches.Circle(
            (cx, cy), r_px, linewidth=0.5, edgecolor="yellow", facecolor="none"
        ))
        # Radius line so you can see exactly what's being measured
        axes[1].plot([cx, cx + r_px], [cy, cy], color="yellow", linewidth=0.4)

    axes[1].set_title(f"Gold NPs detected (n={len(df)})")
    axes[1].axis("off")

    if title:
        fig.suptitle(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_histogram(df_all: pd.DataFrame, out_path: Path) -> None:
    data = df_all["diameter_nm"]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(data, bins=30, edgecolor="white", color="gold")
    ax.axvline(data.mean(), color="tomato", linewidth=2,
               label=f"mean = {data.mean():.1f} ± {data.std():.1f} nm  (n={len(data)})")
    ax.set_xlabel("Gold NP diameter (nm)")
    ax.set_ylabel("Count")
    ax.set_title("Gold NP size distribution")
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
    parser.add_argument("image_path", type=Path,
                        help="Image file (.dm3/.dm4) or folder of images")
    parser.add_argument("--pattern", default="*",
                        help="Glob when image_path is a folder (default: all dm files)")
    parser.add_argument("--out", type=Path, default=Path("results"),
                        help="Output directory (default: results/)")
    parser.add_argument("--gold-threshold", type=float, default=None,
                        help="Intensity cutoff for gold detection (default: auto per image; "
                             "set manually if auto gives wrong results)")
    parser.add_argument("--min-gold-nm", type=float, default=10.0,
                        help="Minimum NP diameter to keep (default 10 nm)")
    parser.add_argument("--max-gold-nm", type=float, default=30.0,
                        help="Maximum NP diameter to keep (default 30 nm)")
    args = parser.parse_args()

    p = args.image_path
    if p.is_file():
        files = [p]
    elif p.is_dir():
        files = sorted(
            f for ext in (".dm3", ".dm4")
            for f in p.glob(f"{args.pattern}{ext}")
        )
    else:
        print(f"Path not found: {p}", file=sys.stderr)
        sys.exit(1)

    if not files:
        print("No .dm3 / .dm4 files found.", file=sys.stderr)
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
        df, t_used = detect_gold(img_norm, nm_per_px, args.gold_threshold,
                                 args.min_gold_nm, args.max_gold_nm)
        df.insert(0, "file", path.name)
        all_rows.append(df)

        n = len(df)
        if n:
            print(f"  threshold  {t_used:.4f} (auto)" if args.gold_threshold is None
                  else f"  threshold  {t_used:.4f} (manual)")
            print(f"  {n} gold NPs  |  mean {df['diameter_nm'].mean():.1f} ± "
                  f"{df['diameter_nm'].std():.1f} nm")
        else:
            print("  No gold NPs detected — try adjusting --gold-threshold or --min/max-gold-nm")

        overlay_path = overlay_dir / (path.stem + "_overlay.png")
        save_overlay(img_norm, df, nm_per_px, overlay_path, title=path.name)
        print(f"  Overlay    → {overlay_path}")

    results = pd.concat(all_rows, ignore_index=True)

    csv_path = args.out / "gold_np_diameters.csv"
    results.to_csv(csv_path, index=False)
    print(f"\n  CSV        → {csv_path}")

    hist_path = args.out / "gold_np_histogram.png"
    save_histogram(results, hist_path)

    print("\n── Summary ──────────────────────────────────────────────")
    data = results["diameter_nm"]
    print(f"  n       {len(data)}")
    print(f"  mean    {data.mean():.1f} nm")
    print(f"  std     {data.std():.1f} nm")
    print(f"  median  {data.median():.1f} nm")
    print(f"  range   {data.min():.1f} – {data.max():.1f} nm")


if __name__ == "__main__":
    main()
