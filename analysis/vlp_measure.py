"""
vlp_measure.py — Gold NP + capsid measurement in VLP TEM images

Detects gold nanoparticles (dense, near-black circular blobs) in Gatan
.dm3 / .dm4 TEM images, then finds the surrounding capsid shell via
radial intensity profile gradient detection.

Two measurements per particle:
  - gold_diameter_nm   : diameter of the gold nanoparticle core
  - capsid_diameter_nm : outer diameter of the capsid shell
                         (NaN for naked cores with no detectable capsid)

Usage:
    uv run python vlp_measure.py <image_path> [options]

Examples:
    uv run python vlp_measure.py "images/VLPs for machine learning project/VLP17_0001.dm4"
    uv run python vlp_measure.py "images/VLPs for machine learning project" --pattern "VLP17_*"
    uv run python vlp_measure.py "images/" --debug
    uv run python vlp_measure.py "images/" --max-capsid-nm 80
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
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
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
    Uses Otsu on the dark half of the histogram. Falls back to manual value
    if provided.
    """
    if gold_threshold is not None:
        return gold_threshold

    dark_pixels = img_norm[img_norm < 0.5]
    if len(dark_pixels) < 100:
        return 0.05

    try:
        return float(filters.threshold_otsu(dark_pixels))
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
    Returns (DataFrame: centroid_y, centroid_x, gold_diameter_nm, gold_radius_px),
    threshold_used.
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
        return pd.DataFrame(
            columns=["centroid_y", "centroid_x", "gold_diameter_nm", "gold_radius_px"]
        ), t

    df["gold_diameter_nm"] = 2 * np.sqrt(df["area"] / np.pi) * nm_per_px
    df["gold_radius_px"] = np.sqrt(df["area"] / np.pi)

    df = df[
        (df["solidity"] > 0.75)
        & (df["eccentricity"] < 0.75)
        & (df["gold_diameter_nm"] >= min_gold_nm)
        & (df["gold_diameter_nm"] <= max_gold_nm)
    ]

    df = df.rename(columns={"centroid-0": "centroid_y", "centroid-1": "centroid_x"})
    return (
        df[["centroid_y", "centroid_x", "gold_diameter_nm", "gold_radius_px"]].reset_index(drop=True),
        t,
    )


# ── Capsid detection via radial profile ───────────────────────────────────────

def radial_profile(
    img: np.ndarray,
    cy: float,
    cx: float,
    r_start: int,
    r_end: int,
    n_angles: int = 360,
) -> np.ndarray:
    """
    Mean intensity at each radius from r_start to r_end pixels from (cy, cx).
    Returns array of length (r_end - r_start).
    """
    h, w = img.shape
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    profile = np.zeros(r_end - r_start)

    for i, r in enumerate(range(r_start, r_end)):
        ys = np.clip((cy + r * np.sin(angles)).astype(int), 0, h - 1)
        xs = np.clip((cx + r * np.cos(angles)).astype(int), 0, w - 1)
        profile[i] = img[ys, xs].mean()

    return profile


def find_capsid_edge(
    img_norm: np.ndarray,
    cy: float,
    cx: float,
    gold_radius_px: float,
    max_capsid_nm: float,
    nm_per_px: float,
    smooth_sigma: float = 3.0,
) -> float | None:
    """
    Find capsid outer radius via radial intensity profile gradient.

    Searches outward from just beyond the gold NP edge. The capsid outer wall
    is where the intensity rises most steeply (peak gradient) as it transitions
    from the darker capsid region to the bright background.

    smooth_sigma controls Gaussian smoothing of the profile before gradient
    computation — higher values are more robust to the speckled capsid texture
    but may shift the detected edge slightly.

    Returns capsid radius in pixels, or None if no clear edge found.
    """
    r_start = max(1, int(gold_radius_px))
    r_end = r_start + int(max_capsid_nm / nm_per_px)
    r_end = min(r_end, min(img_norm.shape) // 2)

    if r_end <= r_start:
        return None

    profile = radial_profile(img_norm, cy, cx, r_start, r_end)

    # Smooth to suppress speckle noise before finding the minimum
    smoothed = gaussian_filter1d(profile, sigma=smooth_sigma)

    # The capsid wall is a dark ring = local minimum in the radial profile.
    # Skip the first few nm beyond the gold edge to avoid catching the
    # intensity recovery from the gold NP itself as a false minimum.
    skip_nm = 5.0
    skip_px = int(skip_nm / nm_per_px)

    search = smoothed[skip_px:]
    if len(search) < 3:
        return None

    # Find minima as peaks in the inverted profile.
    # prominence: the dip must stand out clearly from its surroundings.
    # distance: minima must be at least 3 nm apart (avoids noise spikes).
    min_dist_px = max(1, int(3.0 / nm_per_px))
    peaks, props = find_peaks(-search, prominence=0.02, distance=min_dist_px)

    if len(peaks) == 0:
        return None

    peak_i = peaks[0]

    # The dark stain ring minimum is just outside the capsid protein.
    # The true outer capsid boundary is the steepest descending gradient
    # on the approach to the minimum (where bright protein → dark stain).
    approach = search[: peak_i + 1]
    if len(approach) >= 3:
        grad = np.gradient(approach)
        edge_i = int(np.argmin(grad))  # most negative slope = outer protein wall
        # Sub-pixel refinement on the gradient minimum
        if 0 < edge_i < len(grad) - 1:
            g0, g1, g2 = grad[edge_i - 1], grad[edge_i], grad[edge_i + 1]
            denom = 2 * (g0 - 2 * g1 + g2)
            if denom != 0:
                subpixel_offset = (g0 - g2) / denom
                return float(r_start + skip_px + edge_i + subpixel_offset)
        return float(r_start + skip_px + edge_i)

    # Fallback: return the minimum itself
    if 0 < peak_i < len(search) - 1:
        y0, y1, y2 = search[peak_i - 1], search[peak_i], search[peak_i + 1]
        denom = 2 * (y0 - 2 * y1 + y2)
        if denom != 0:
            subpixel_offset = (y0 - y2) / denom
            return float(r_start + skip_px + peak_i + subpixel_offset)

    return float(r_start + skip_px + peak_i)


# ── Output ────────────────────────────────────────────────────────────────────

def save_overlay(
    img_norm: np.ndarray,
    df: pd.DataFrame,
    nm_per_px: float,
    out_path: Path,
    title: str = "",
    highlight_indices: set | None = None,
) -> None:
    """
    Gold NP in yellow, capsid in cyan, radius lines for both.
    Particles in highlight_indices are drawn in orange so you can match
    them to the debug profile plots.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    axes[0].imshow(img_norm, cmap="gray", interpolation="none")
    axes[0].set_title("Image")
    axes[0].axis("off")

    axes[1].imshow(img_norm, cmap="gray", interpolation="none")
    for idx, row in df.iterrows():
        cx, cy = row["centroid_x"], row["centroid_y"]
        highlighted = highlight_indices is not None and idx in highlight_indices
        gold_color = "orange" if highlighted else "yellow"

        # Gold NP
        r_gold = row["gold_diameter_nm"] / 2 / nm_per_px
        axes[1].add_patch(mpatches.Circle(
            (cx, cy), r_gold, linewidth=0.8 if highlighted else 0.5,
            edgecolor=gold_color, facecolor="none"
        ))
        axes[1].plot([cx, cx + r_gold], [cy, cy], color=gold_color, linewidth=0.4)

        # Capsid (only if detected)
        if not np.isnan(row["capsid_diameter_nm"]):
            r_cap = row["capsid_diameter_nm"] / 2 / nm_per_px
            axes[1].add_patch(mpatches.Circle(
                (cx, cy), r_cap, linewidth=0.8 if highlighted else 0.5,
                edgecolor="orange" if highlighted else "cyan", facecolor="none"
            ))
            axes[1].plot([cx, cx + r_cap], [cy, cy],
                         color="orange" if highlighted else "cyan", linewidth=0.4)

        # Index label — always shown, larger for highlighted particles
        axes[1].text(
            cx, cy - r_gold - 2, str(idx),
            color="orange" if highlighted else "white",
            fontsize=5 if highlighted else 3.5,
            fontweight="bold" if highlighted else "normal",
            ha="center", va="bottom",
            bbox=dict(boxstyle="round,pad=0.1", facecolor="black", alpha=0.6, linewidth=0),
        )

    n_total = len(df)
    n_capsid = df["capsid_diameter_nm"].notna().sum()
    axes[1].set_title(f"n={n_total}  ({n_capsid} with capsid)")
    axes[1].axis("off")
    axes[1].legend(handles=[
        mpatches.Patch(color="yellow", label="Gold NP"),
        mpatches.Patch(color="cyan", label="Capsid"),
    ], loc="lower right", fontsize=7)

    if title:
        fig.suptitle(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_debug_profiles(
    img_norm: np.ndarray,
    df: pd.DataFrame,
    nm_per_px: float,
    max_capsid_nm: float,
    out_path: Path,
    indices: list[int] | None = None,
) -> None:
    """Plot radial profiles for specific particles so you can inspect capsid detection."""
    if indices is not None:
        sample = df.loc[[i for i in indices if i in df.index]]
    else:
        sample = df.dropna(subset=["capsid_diameter_nm"]).head(6)
    if sample.empty:
        return

    cols = min(len(sample), 3)
    rows = (len(sample) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
    axes = np.array(axes).flatten()

    for ax, (idx, row) in zip(axes, sample.iterrows()):
        cy, cx = row["centroid_y"], row["centroid_x"]
        gold_r = row["gold_radius_px"]
        r_start = max(1, int(gold_r))
        r_end = r_start + int(max_capsid_nm / nm_per_px)
        r_end = min(r_end, min(img_norm.shape) // 2)

        profile = radial_profile(img_norm, cy, cx, r_start, r_end)
        smoothed = gaussian_filter1d(profile, sigma=3.0)
        radii_nm = (r_start + np.arange(len(profile))) * nm_per_px

        ax.plot(radii_nm, profile, color="lightsteelblue", linewidth=1, label="raw")
        ax.plot(radii_nm, smoothed, color="steelblue", linewidth=1.5, label="smoothed")
        ax.axvline(row["gold_diameter_nm"] / 2, color="yellow", linewidth=1.5,
                   label=f"gold r = {row['gold_diameter_nm']/2:.1f} nm")
        ax.axvline(row["capsid_diameter_nm"] / 2, color="cyan", linewidth=1.5,
                   label=f"capsid r = {row['capsid_diameter_nm']/2:.1f} nm")
        ax.set_xlabel("Radius (nm)")
        ax.set_ylabel("Intensity")
        ax.set_title(f"#{idx}  ({cx:.0f}, {cy:.0f})")
        ax.legend(fontsize=7)

    for ax in axes[len(sample):]:
        ax.set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Debug      → {out_path}")


def save_scatter(df_all: pd.DataFrame, out_path: Path) -> None:
    df = df_all.dropna(subset=["capsid_diameter_nm"])
    reliable = df[df["is_reliable"]] if "is_reliable" in df.columns else df
    unreliable = df[~df["is_reliable"]] if "is_reliable" in df.columns else pd.DataFrame()

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(reliable["gold_diameter_nm"], reliable["capsid_diameter_nm"],
               alpha=0.5, s=18, color="steelblue", edgecolors="none", label="reliable")
    if not unreliable.empty:
        ax.scatter(unreliable["gold_diameter_nm"], unreliable["capsid_diameter_nm"],
                   alpha=0.7, s=40, color="tomato", marker="x",
                   label=f"flagged ({len(unreliable)})")
    ax.axvline(reliable["gold_diameter_nm"].mean(), color="gold", linewidth=1.2, linestyle="--",
               label=f"gold mean = {reliable['gold_diameter_nm'].mean():.1f} nm")
    ax.axhline(reliable["capsid_diameter_nm"].mean(), color="cyan", linewidth=1.2, linestyle="--",
               label=f"capsid mean = {reliable['capsid_diameter_nm'].mean():.1f} nm")
    ax.set_xlabel("Gold NP diameter (nm)")
    ax.set_ylabel("VLP capsid diameter (nm)")
    ax.set_title(f"VLP capsid vs gold NP diameter  (n={len(reliable)} reliable, {len(unreliable)} flagged)")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Scatter    → {out_path}")


def save_histograms(df_all: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    for ax, col, label, color in [
        (axes[0], "gold_diameter_nm", "Gold NP diameter (nm)", "gold"),
        (axes[1], "capsid_diameter_nm", "Capsid diameter (nm)", "steelblue"),
    ]:
        data = df_all[col].dropna()
        if len(data) == 0:
            ax.set_title(f"{label}\n(no data)")
            continue
        bins = np.arange(np.floor(data.min()), np.ceil(data.max()) + 1, 1)
        ax.hist(data, bins=bins, edgecolor="white", linewidth=0.5, color=color)
        ax.axvline(data.mean(), color="tomato", linewidth=2,
                   label=f"mean = {data.mean():.1f} ± {data.std():.1f} nm\n(n={len(data)})")
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
                        help="Glob when image_path is a folder (default: all dm files)")
    parser.add_argument("--out", type=Path, default=Path("results"),
                        help="Output directory (default: results/)")

    # Gold NP tuning
    parser.add_argument("--gold-threshold", type=float, default=None,
                        help="Intensity cutoff for gold detection (default: auto per image)")
    parser.add_argument("--min-gold-nm", type=float, default=10.0,
                        help="Minimum gold NP diameter (default 10 nm)")
    parser.add_argument("--max-gold-nm", type=float, default=30.0,
                        help="Maximum gold NP diameter (default 30 nm)")

    # Capsid tuning
    parser.add_argument("--max-capsid-nm", type=float, default=150.0,
                        help="Max distance from gold edge to search for capsid (default 150 nm)")
    parser.add_argument("--smooth-sigma", type=float, default=3.0,
                        help="Gaussian smoothing of radial profile before gradient (default 3.0; "
                             "increase if capsid ring is very noisy)")

    # Debug
    parser.add_argument("--debug", action="store_true",
                        help="Save radial profile plots to inspect capsid edge detection")
    parser.add_argument("--debug-indices", type=int, nargs="+", default=None,
                        help="Specific particle indices to show in debug profiles (default: first 6 with capsid)")

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

        df_gold, t_used = detect_gold(
            img_norm, nm_per_px, args.gold_threshold, args.min_gold_nm, args.max_gold_nm
        )
        print(f"  threshold  {t_used:.4f} {'(auto)' if args.gold_threshold is None else '(manual)'}")
        print(f"  {len(df_gold)} gold NPs detected")

        # Measure capsid for each gold NP
        capsid_diameters = []
        capsid_radii_px = []
        for _, g in df_gold.iterrows():
            r_px = find_capsid_edge(
                img_norm,
                g["centroid_y"], g["centroid_x"],
                g["gold_radius_px"],
                args.max_capsid_nm,
                nm_per_px,
                smooth_sigma=args.smooth_sigma,
            )
            capsid_radii_px.append(r_px)
            capsid_diameters.append(2 * r_px * nm_per_px if r_px else np.nan)

        df_gold["capsid_diameter_nm"] = capsid_diameters
        df_gold["capsid_radius_px"] = capsid_radii_px
        df_gold = df_gold.reset_index(drop=True)  # ensure index matches overlay labels

        # Reliability flag: capsid within 2 std of per-image median
        cd = df_gold["capsid_diameter_nm"].dropna()
        if len(cd) > 2:
            median = cd.median()
            std = cd.std()
            df_gold["is_reliable"] = (
                df_gold["capsid_diameter_nm"].isna() |
                df_gold["capsid_diameter_nm"].between(median - 2 * std, median + 2 * std)
            )
        else:
            df_gold["is_reliable"] = True

        df_gold.insert(0, "file", path.name)
        all_rows.append(df_gold)

        n_capsid = df_gold["capsid_diameter_nm"].notna().sum()
        n_unreliable = (~df_gold["is_reliable"]).sum()
        print(f"  {n_capsid}/{len(df_gold)} capsids detected  |  {n_unreliable} flagged unreliable")
        if n_capsid:
            reliable = df_gold[df_gold["is_reliable"] & df_gold["capsid_diameter_nm"].notna()]
            print(f"  gold   mean {df_gold['gold_diameter_nm'].mean():.1f} ± {df_gold['gold_diameter_nm'].std():.1f} nm")
            print(f"  capsid mean {reliable['capsid_diameter_nm'].mean():.1f} ± {reliable['capsid_diameter_nm'].std():.1f} nm  (reliable only)")

        # Determine which particles will appear in debug profiles
        highlight = None
        if args.debug:
            if args.debug_indices:
                sample_idx = [i for i in args.debug_indices if i in df_gold.index]
            else:
                sample_idx = df_gold.dropna(subset=["capsid_diameter_nm"]).head(6).index.tolist()
            highlight = set(sample_idx)

        overlay_path = overlay_dir / (path.stem + "_overlay.png")
        save_overlay(img_norm, df_gold, nm_per_px, overlay_path,
                     title=path.name, highlight_indices=highlight)
        print(f"  Overlay    → {overlay_path}")

        if args.debug:
            debug_path = overlay_dir / (path.stem + "_profiles.png")
            save_debug_profiles(img_norm, df_gold, nm_per_px, args.max_capsid_nm,
                                debug_path, indices=sample_idx)

    results = pd.concat(all_rows, ignore_index=True)

    csv_path = args.out / "vlp_measurements.csv"
    results.to_csv(csv_path, index=False)
    print(f"\n  CSV        → {csv_path}")

    hist_path = args.out / "vlp_histograms.png"
    save_histograms(results, hist_path)

    scatter_path = args.out / "vlp_scatter.png"
    save_scatter(results, scatter_path)

    print("\n── Summary ──────────────────────────────────────────────")
    gold = results["gold_diameter_nm"].dropna()
    print(f"\nGold NP (n={len(gold)}):")
    print(f"  mean    {gold.mean():.1f} nm")
    print(f"  std     {gold.std():.1f} nm")
    print(f"  median  {gold.median():.1f} nm")
    print(f"  range   {gold.min():.1f} – {gold.max():.1f} nm")

    reliable_capsid = results[results["is_reliable"]]["capsid_diameter_nm"].dropna()
    all_capsid = results["capsid_diameter_nm"].dropna()
    n_flagged = (~results["is_reliable"]).sum()
    print(f"\nCapsid — reliable only (n={len(reliable_capsid)}, {n_flagged} flagged):")
    print(f"  mean    {reliable_capsid.mean():.1f} nm")
    print(f"  std     {reliable_capsid.std():.1f} nm")
    print(f"  median  {reliable_capsid.median():.1f} nm")
    print(f"  range   {reliable_capsid.min():.1f} – {reliable_capsid.max():.1f} nm")


if __name__ == "__main__":
    main()
