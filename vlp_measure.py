"""
vlp_measure.py — VLP size analysis from TEM images

Detects virus-like particles (VLPs) in Gatan .dm3 / .dm4 TEM images.
Each VLP consists of a dense gold nanoparticle core (very dark, solid circle)
surrounded by a speckled capsid ring.

Two measurements per particle:
  - gold_diameter_nm  : diameter of the gold nanoparticle core
  - capsid_diameter_nm: outer diameter of the capsid shell
                        (NaN for naked cores with no detectable capsid)

Strategy:
  1. Detect gold NPs first — highest contrast feature, easy anchor point.
  2. For each gold NP center, extract a radial intensity profile outward.
  3. Find where the profile recovers to background level → capsid outer edge.

Usage:
    uv run python vlp_measure.py <image_dir> [options]

Examples:
    uv run python vlp_measure.py "images/VLPs for machine learning project"
    uv run python vlp_measure.py "images/" --pattern "VLP17_*" --out results/
    uv run python vlp_measure.py "images/VLP17_0001.dm4" --debug
    uv run python vlp_measure.py "images/" --capsid-threshold 0.75
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
from skimage import exposure, morphology, measure


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

def detect_gold(
    img_norm: np.ndarray,
    nm_per_px: float,
    threshold: float,
) -> pd.DataFrame:
    """
    Find gold nanoparticles by intensity thresholding.

    Gold NPs are the darkest feature in the image (near-black on TEM).
    threshold is the intensity cutoff in [0,1] — pixels below this are
    considered gold. Default 0.05 means bottom 5% of the normalised range.
    Lower if gold NPs are being missed; raise if background noise is detected.

    Returns DataFrame with columns: centroid_y, centroid_x, gold_radius_px,
    gold_diameter_nm.
    """
    # Threshold: keep only the darkest pixels
    binary = img_norm < threshold

    # Clean up single-pixel noise
    binary = morphology.remove_small_objects(binary, min_size=4)

    # Label connected regions
    labels = ndi.label(binary)[0]

    props = measure.regionprops_table(
        labels,
        properties=("label", "area", "centroid", "eccentricity", "solidity"),
    )
    df = pd.DataFrame(props)
    if df.empty:
        return pd.DataFrame(
            columns=["centroid_y", "centroid_x", "gold_radius_px", "gold_diameter_nm"]
        )

    # Equivalent circular diameter from area
    df["gold_radius_px"] = np.sqrt(df["area"] / np.pi)
    df["gold_diameter_nm"] = 2 * df["gold_radius_px"] * nm_per_px

    # Keep only reasonably circular regions (gold NPs are round)
    df = df[(df["solidity"] > 0.75) & (df["eccentricity"] < 0.75)]

    df = df.rename(columns={"centroid-0": "centroid_y", "centroid-1": "centroid_x"})
    return df[["centroid_y", "centroid_x", "gold_radius_px", "gold_diameter_nm"]].reset_index(drop=True)


# ── Radial profile + capsid measurement ───────────────────────────────────────

def radial_profile(
    img: np.ndarray,
    cy: float,
    cx: float,
    max_radius_px: int,
    n_angles: int = 360,
) -> np.ndarray:
    """
    Sample image intensity at evenly spaced radii from (cy, cx) outward,
    averaging over n_angles angular samples at each radius.
    Returns 1-D array of mean intensity vs radius (length = max_radius_px).
    """
    radii = np.arange(1, max_radius_px + 1)
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    profile = np.zeros(max_radius_px)

    h, w = img.shape
    for i, r in enumerate(radii):
        ys = (cy + r * np.sin(angles)).astype(int)
        xs = (cx + r * np.cos(angles)).astype(int)
        # Keep only samples inside image bounds
        valid = (ys >= 0) & (ys < h) & (xs >= 0) & (xs < w)
        if valid.sum() > 0:
            profile[i] = img[ys[valid], xs[valid]].mean()

    return profile


def find_capsid_radius(
    profile: np.ndarray,
    gold_radius_px: float,
    capsid_threshold: float,
) -> float | None:
    """
    Find the capsid outer radius from a radial intensity profile.

    The profile starts dark (gold NP), rises through the capsid, then levels
    off at background. The capsid outer edge is where the profile first crosses
    capsid_threshold * background_level on its way up.

    Returns radius in pixels, or None if no clear capsid edge is found.
    """
    # Estimate background as the mean of the outermost 20% of the profile
    tail_start = int(len(profile) * 0.80)
    background = profile[tail_start:].mean()

    # Target crossing level
    target = capsid_threshold * background

    # Search outward from just beyond the gold NP
    start = max(0, int(gold_radius_px))
    for i in range(start, len(profile)):
        if profile[i] >= target:
            return float(i)

    return None  # no clear capsid edge found


# ── Per-file processing ───────────────────────────────────────────────────────

def process_image(
    path: Path,
    gold_threshold: float,
    capsid_threshold: float,
    max_search_nm: float,
    debug: bool,
    debug_n: int,
) -> pd.DataFrame:
    """Detect VLPs in one image and return a DataFrame of measurements."""
    img, nm_per_px = load_dm(path)
    print(f"  {img.shape[0]}×{img.shape[1]} px  |  {nm_per_px:.4f} nm/px")

    img_norm = normalise(img)

    # 1. Detect gold NPs
    df_gold = detect_gold(img_norm, nm_per_px, threshold=gold_threshold)
    print(f"  {len(df_gold)} gold NPs detected")

    if df_gold.empty:
        print("  → Try lowering --gold-threshold")
        return pd.DataFrame()

    # 2. For each gold NP, measure capsid via radial profile
    max_search_px = int(max_search_nm / nm_per_px)
    rows = []
    debug_particles = []

    for _, g in df_gold.iterrows():
        cy, cx = g["centroid_y"], g["centroid_x"]
        gold_r = g["gold_radius_px"]

        profile = radial_profile(img_norm, cy, cx, max_search_px)
        capsid_r = find_capsid_radius(profile, gold_r, capsid_threshold)

        row = {
            "file": path.name,
            "centroid_y": cy,
            "centroid_x": cx,
            "gold_diameter_nm": g["gold_diameter_nm"],
            "gold_radius_px": gold_r,
            "capsid_diameter_nm": 2 * capsid_r * nm_per_px if capsid_r else np.nan,
            "capsid_radius_px": capsid_r if capsid_r else np.nan,
            "has_capsid": capsid_r is not None,
        }
        rows.append(row)

        if debug and len(debug_particles) < debug_n:
            debug_particles.append((profile, gold_r, capsid_r, cy, cx, g["gold_diameter_nm"]))

    df = pd.DataFrame(rows)

    n_capsid = df["has_capsid"].sum()
    n_naked = (~df["has_capsid"]).sum()
    print(f"  {n_capsid} with capsid  |  {n_naked} naked cores")
    if n_capsid > 0:
        print(f"  Gold mean:   {df['gold_diameter_nm'].mean():.1f} ± {df['gold_diameter_nm'].std():.1f} nm")
        print(f"  Capsid mean: {df.loc[df['has_capsid'], 'capsid_diameter_nm'].mean():.1f} ± "
              f"{df.loc[df['has_capsid'], 'capsid_diameter_nm'].std():.1f} nm")

    return df, img_norm, nm_per_px, debug_particles


# ── Output ────────────────────────────────────────────────────────────────────

def save_overlay(
    img_norm: np.ndarray,
    df: pd.DataFrame,
    nm_per_px: float,
    out_path: Path,
    title: str = "",
) -> None:
    """Draw two concentric circles per particle: gold (yellow) and capsid (cyan)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    axes[0].imshow(img_norm, cmap="gray", interpolation="none")
    axes[0].set_title("Image")
    axes[0].axis("off")

    axes[1].imshow(img_norm, cmap="gray", interpolation="none")
    for _, row in df.iterrows():
        cy, cx = row["centroid_y"], row["centroid_x"]

        # Gold NP circle (yellow)
        r_gold = row["gold_radius_px"]
        axes[1].add_patch(mpatches.Circle(
            (cx, cy), r_gold, linewidth=0.8, edgecolor="yellow", facecolor="none"
        ))
        # Radius line so you can see exactly what's being measured
        axes[1].plot([cx, cx + r_gold], [cy, cy], color="yellow", linewidth=0.6)

        # Capsid circle (cyan) — only if detected
        if row["has_capsid"] and not np.isnan(row["capsid_radius_px"]):
            r_cap = row["capsid_radius_px"]
            axes[1].add_patch(mpatches.Circle(
                (cx, cy), r_cap, linewidth=0.8, edgecolor="cyan", facecolor="none"
            ))
            axes[1].plot([cx, cx + r_cap], [cy, cy], color="cyan", linewidth=0.6)

    n_total = len(df)
    n_capsid = df["has_capsid"].sum()
    axes[1].set_title(f"n={n_total} particles  ({n_capsid} with capsid)")
    axes[1].axis("off")

    legend = [
        mpatches.Patch(color="yellow", label="Gold NP"),
        mpatches.Patch(color="cyan", label="Capsid"),
    ]
    axes[1].legend(handles=legend, loc="lower right", fontsize=8)

    if title:
        fig.suptitle(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_debug_profiles(debug_particles: list, nm_per_px: float, out_path: Path) -> None:
    """Plot radial intensity profiles for a handful of particles."""
    n = len(debug_particles)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]

    for ax, (profile, gold_r, capsid_r, cy, cx, gold_d) in zip(axes, debug_particles):
        radii_nm = np.arange(len(profile)) * nm_per_px
        ax.plot(radii_nm, profile, color="steelblue", linewidth=1.5)
        ax.axvline(gold_r * nm_per_px, color="yellow", linewidth=1.5,
                   label=f"gold r = {gold_r * nm_per_px:.1f} nm")
        if capsid_r:
            ax.axvline(capsid_r * nm_per_px, color="cyan", linewidth=1.5,
                       label=f"capsid r = {capsid_r * nm_per_px:.1f} nm")
        ax.set_xlabel("Radius (nm)")
        ax.set_ylabel("Intensity")
        ax.set_title(f"Particle at ({cx:.0f}, {cy:.0f})")
        ax.legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Debug profiles → {out_path}")


def save_histograms(df_all: pd.DataFrame, out_path: Path) -> None:
    """Side-by-side histograms: gold diameter and capsid diameter."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    for ax, col, label, color in [
        (axes[0], "gold_diameter_nm", "Gold NP diameter (nm)", "gold"),
        (axes[1], "capsid_diameter_nm", "Capsid diameter (nm)", "steelblue"),
    ]:
        data = df_all[col].dropna()
        if len(data) == 0:
            ax.set_title(f"{label}\n(no data)")
            continue
        ax.hist(data, bins=30, edgecolor="white", color=color)
        mean_d = data.mean()
        std_d = data.std()
        ax.axvline(mean_d, color="tomato", linewidth=2,
                   label=f"mean = {mean_d:.1f} ± {std_d:.1f} nm\n(n={len(data)})")
        ax.set_xlabel(label)
        ax.set_ylabel("Count")
        ax.legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Histograms → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("image_path", type=Path,
                        help="Image file (.dm3/.dm4) or folder of images")
    parser.add_argument("--pattern", default="*",
                        help="Glob pattern when image_path is a folder (default: all dm files)")
    parser.add_argument("--out", type=Path, default=Path("results"),
                        help="Output directory (default: results/)")

    # Detection tuning
    parser.add_argument("--gold-threshold", type=float, default=0.05,
                        help="Gold NP detection sensitivity (default 0.05; lower = more detections)")
    parser.add_argument("--capsid-threshold", type=float, default=0.70,
                        help="Fraction of background intensity that marks the capsid edge "
                             "(default 0.70; raise if capsid boundary is being set too far out)")
    parser.add_argument("--max-search-nm", type=float, default=300,
                        help="How far outward to search for the capsid edge (default 300 nm)")

    # Debug
    parser.add_argument("--debug", action="store_true",
                        help="Save radial profile plots so you can inspect what's being measured")
    parser.add_argument("--debug-n", type=int, default=4,
                        help="Number of particles to show in debug plots (default 4)")

    args = parser.parse_args()

    # Collect files
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
        print(f"No .dm3 / .dm4 files found.", file=sys.stderr)
        sys.exit(1)

    args.out.mkdir(parents=True, exist_ok=True)
    overlay_dir = args.out / "overlays"
    overlay_dir.mkdir(exist_ok=True)

    all_rows = []

    for path in files:
        print(f"\n{path.name}")
        result = process_image(
            path,
            gold_threshold=args.gold_threshold,
            capsid_threshold=args.capsid_threshold,
            max_search_nm=args.max_search_nm,
            debug=args.debug,
            debug_n=args.debug_n,
        )
        if isinstance(result, tuple):
            df, img_norm, nm_per_px, debug_particles = result
        else:
            continue

        if df.empty:
            continue

        all_rows.append(df)

        overlay_path = overlay_dir / (path.stem + "_overlay.png")
        save_overlay(img_norm, df, nm_per_px, overlay_path, title=path.name)
        print(f"  Overlay    → {overlay_path}")

        if args.debug and debug_particles:
            debug_path = overlay_dir / (path.stem + "_debug_profiles.png")
            save_debug_profiles(debug_particles, nm_per_px, debug_path)

    if not all_rows:
        print("\nNo particles detected.")
        return

    results = pd.concat(all_rows, ignore_index=True)

    csv_path = args.out / "vlp_measurements.csv"
    results.to_csv(csv_path, index=False)
    print(f"\n  CSV        → {csv_path}")

    hist_path = args.out / "vlp_histograms.png"
    save_histograms(results, hist_path)

    print("\n── Summary ──────────────────────────────────────────────")
    for col, label in [("gold_diameter_nm", "Gold NP"), ("capsid_diameter_nm", "Capsid")]:
        data = results[col].dropna()
        if len(data):
            print(f"\n{label} (n={len(data)}):")
            print(f"  mean   {data.mean():.1f} nm")
            print(f"  std    {data.std():.1f} nm")
            print(f"  median {data.median():.1f} nm")
            print(f"  range  {data.min():.1f} – {data.max():.1f} nm")


if __name__ == "__main__":
    main()
