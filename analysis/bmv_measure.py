"""
bmv_measure.py — BMV / BOG capsid measurement in TEM images

Negative-stained BMV (and fluorescent-protein-labeled BMV "BOG") capsids
appear as donut-shaped objects: a dark stain-filled center, a bright
protein-wall ring, and a dark uranyl-stain pool around the outside. Unlike
gold-NP VLPs there is no high-contrast anchor at the center, so this script
detects whole capsids by their circular-edge signature.

Pipeline:
  1. Hough circle transform at the expected ~28 nm radius (works for
     donut signature regardless of bright/dark interior).
  2. Center refinement: snap each candidate to the local intensity minimum
     (dark stain core) — but only if the snap improves wall circularity.
  3. Wall fit: radial intensity profile from each center; the wall is
     placed where intensity has fallen 75% of the way from the bright-ring
     peak to the post-peak stain minimum (≈ visible outer protein edge).
  4. Quality filtering: bright-ring vs dark-stain contrast, per-sector
     wall-radius circularity, absolute exterior darkness.
  5. Overlap exclusion: drop reliable pairs whose centres are closer than
     (r1 + r2) — touching particles share a stain ring that contaminates
     both wall fits.

All tuning constants live as named globals at the top of this file. The
profile smoothing constant is in nm and converted to pixels per-image so
detections are magnification-invariant.

Usage:
    uv run python bmv_measure.py images/BMV/BMV_a.dm3
    uv run python bmv_measure.py images/BMV --pattern "BMV_*" --workers 6
    uv run python bmv_measure.py images/BMV/BMV_a.dm3 --debug
"""

import argparse
import queue
import sys
import threading
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import ncempy.io as nio
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.spatial import cKDTree
from skimage import exposure, filters
from skimage.transform import hough_circle, hough_circle_peaks


# ── Shared utilities (image loading, normalisation, radial profile sampling) ──
# Inlined here so this script is self-contained — drop just this file
# alongside images and run it. Same definitions as in vlp_measure_v2.py.

def load_dm(path: Path) -> tuple[np.ndarray, float]:
    """Return (image as float32, nm_per_pixel) for a Gatan .dm3 / .dm4 file."""
    d = nio.read(str(path))
    img = d["data"].astype(np.float32)
    nm_per_px = float(d["pixelSize"][0])
    return img, nm_per_px


def normalise(img: np.ndarray) -> np.ndarray:
    """Percentile-stretch to [0,1] then CLAHE for local contrast."""
    lo, hi = np.percentile(img, [0.5, 99.5])
    img = np.clip((img - lo) / (hi - lo), 0, 1)
    return exposure.equalize_adapthist(img, clip_limit=0.02)


def radial_profile(
    img: np.ndarray,
    cy: float,
    cx: float,
    r_start: int,
    r_end: int,
    n_angles: int = 360,
) -> np.ndarray:
    """Mean intensity at each integer radius from r_start to r_end (px) around (cy, cx)."""
    h, w = img.shape
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    profile = np.zeros(r_end - r_start)
    for i, r in enumerate(range(r_start, r_end)):
        ys = np.clip((cy + r * np.sin(angles)).astype(int), 0, h - 1)
        xs = np.clip((cx + r * np.cos(angles)).astype(int), 0, w - 1)
        profile[i] = img[ys, xs].mean()
    return profile


# ── Tuning constants ─────────────────────────────────────────────────────────
# Edit these to retune detection / wall fitting / quality filtering.
# Defaults were calibrated against user-labeled good (#349, #124) vs
# bad (#48, #153, #574, #313, #347) examples on BMV_a.dm3.

# Pre-smoothing applied ONLY for detection / center refinement (not wall fit).
# Suppresses sub-capsid features (capsomere subunits, RNA texture) that
# become resolved at high magnification and confuse circle detection. The
# value is set just above the typical capsomere scale (~1.5 nm) so the
# whole capsid still reads as a single bright object.
DETECTION_PRESMOOTH_NM = 2.0

# Hough circle detection
HOUGH_DIAMETER_TOL   = 0.20   # ± fraction of expected radius for vote search
HOUGH_THRESHOLD      = 0.25   # vote cutoff as fraction of max accumulator
HOUGH_MIN_DIST_FRAC  = 1.4    # NMS min center separation = capsids touch only

# Radial profile / wall fitting
PROFILE_R_START_FRAC = 0.5    # inner search start (× expected radius)
PROFILE_R_END_FRAC   = 1.6    # outer search extent (× expected radius)
PROFILE_SMOOTH_NM    = 0.7    # nm of smoothing on the radial profile.
                              # Converted to px per-image (sigma = nm / nm_per_px)
                              # so smoothing is mag-invariant. Was 3.0 px fixed,
                              # which biased measurements small at low mag and
                              # large at high mag (see BMV_h–k +1 nm shift).
WALL_DESCENT_FRAC    = 0.75   # 0.5 = half-max (inner protein-ring edge,
                              # too small visually); 1.0 = post-peak min
                              # (outer stain-ring edge, too big visually);
                              # 0.75 = between, matches what user would draw.

# Center refinement (snap to local intensity minimum if it lowers wall_cv)
CENTER_WINDOW_FRAC   = 0.2    # search window size as × expected radius
CENTER_SMOOTH_FRAC   = 0.25   # smoothing sigma as × expected radius

# Quality filters (must all pass for is_reliable=True)
DIAM_TOL_FRAC        = 0.30   # accepted diameter range = expected × (1 ± this)
MIN_CONTRAST         = 0.32   # ring_mean − exterior_mean (dark stain pool)
MIN_UNIFORMITY       = 0.0    # 1 − ring_std/ring_mean (disabled by default)
MAX_WALL_CV          = 0.18   # per-sector wall-radius coefficient of variation
MAX_EXTERIOR_MEAN    = 0.32   # absolute exterior intensity (must be dark stain)
OVERLAP_TOL          = 0.95   # touching-pair exclusion: dist < (r1+r2)·this


# ── Center detection ─────────────────────────────────────────────────────────

def flatten_background(img_norm: np.ndarray, sigma_px: float) -> np.ndarray:
    """Subtract a heavily-blurred copy to remove uneven stain background."""
    bg = gaussian_filter(img_norm, sigma=sigma_px)
    flat = img_norm - bg
    # rescale to [0,1] for downstream blob detection
    lo, hi = np.percentile(flat, [0.5, 99.5])
    return np.clip((flat - lo) / (hi - lo + 1e-9), 0, 1)


def detect_capsid_centers(
    img_norm: np.ndarray,
    nm_per_px: float,
    expected_diameter_nm: float,
    diameter_tol: float = HOUGH_DIAMETER_TOL,
    n_radii: int = 9,
    min_distance_frac: float = HOUGH_MIN_DIST_FRAC,
    peak_threshold: float = HOUGH_THRESHOLD,
    max_candidates: int = 5000,
) -> pd.DataFrame:
    """
    Find capsid centers using the Hough circle transform.

    Hough circles vote from intensity-gradient pixels toward the circle's
    center, so it locates circular structures regardless of whether the
    interior is bright, dark, or donut-shaped — exactly what we need for
    BMV (dark center + bright protein ring + dark stain outside).

    Pipeline:
      1. Canny edge detection on the normalised image.
      2. Hough vote at a small range of radii bracketing expected.
      3. Non-max suppression at min_distance_frac × expected radius.

    diameter_tol: ±tol fractional radius range to vote over.
    peak_threshold: relative vote threshold (fraction of max accumulator).
        Lower → more candidates, more false positives.
    """
    expected_r_px = expected_diameter_nm / 2 / nm_per_px
    r_min = max(2, int(expected_r_px * (1 - diameter_tol)))
    r_max = max(r_min + 1, int(expected_r_px * (1 + diameter_tol)) + 1)
    radii = np.linspace(r_min, r_max, n_radii).astype(int)
    radii = np.unique(radii)

    # Canny edges. sigma scaled to capsid feature scale to avoid noise edges.
    edges = filters.farid(img_norm)  # gradient magnitude
    edge_thresh = np.percentile(edges, 80)
    edge_map = edges > edge_thresh

    accums = hough_circle(edge_map, radii)

    # Pick peaks from vote accumulator. min_xdistance/min_ydistance enforce NMS.
    min_dist = max(1, int(expected_r_px * min_distance_frac))
    accepted_y, accepted_x, accepted_r = [], [], []
    accepted_score = []
    _, cx, cy, r = hough_circle_peaks(
        accums, radii,
        min_xdistance=min_dist, min_ydistance=min_dist,
        threshold=peak_threshold * accums.max(),
        total_num_peaks=max_candidates,
    )
    accepted_y = list(cy)
    accepted_x = list(cx)
    accepted_r = list(r)
    # vote score from accumulator at each peak
    for y, x, rr in zip(cy, cx, r):
        ri = int(np.searchsorted(radii, rr))
        ri = min(ri, len(radii) - 1)
        accepted_score.append(float(accums[ri, y, x]))

    if not accepted_y:
        return pd.DataFrame(
            columns=["centroid_y", "centroid_x", "blob_radius_px", "hough_score"]
        )

    df = pd.DataFrame({
        "centroid_y": np.array(accepted_y, dtype=float),
        "centroid_x": np.array(accepted_x, dtype=float),
        "blob_radius_px": np.array(accepted_r, dtype=float),
        "hough_score": np.array(accepted_score, dtype=float),
    })
    return df.reset_index(drop=True)


# ── Center refinement ────────────────────────────────────────────────────────

def refine_center(
    img_norm: np.ndarray,
    cy: float, cx: float,
    expected_r_px: float,
    window_frac: float = CENTER_WINDOW_FRAC,
    smooth_frac: float = CENTER_SMOOTH_FRAC,
) -> tuple[float, float]:
    """
    Snap an approximate center to the local intensity minimum nearby.

    BMV particles are donut-shaped with a dark center (stain-filled
    interior). Hough vote peaks can sit a few pixels off the true center,
    which biases the wall fit small. We smooth at ~quarter of capsid
    radius and pick the argmin in a tight window so the snap can't escape
    into a neighbouring particle or stain pool.
    """
    h, w = img_norm.shape
    half = max(2, int(expected_r_px * window_frac))
    y0, y1 = max(0, int(cy) - half), min(h, int(cy) + half + 1)
    x0, x1 = max(0, int(cx) - half), min(w, int(cx) + half + 1)
    crop = img_norm[y0:y1, x0:x1]
    if crop.size == 0:
        return cy, cx

    sigma = max(1.0, expected_r_px * smooth_frac)
    crop_smooth = gaussian_filter(crop, sigma=sigma)
    iy, ix = np.unravel_index(crop_smooth.argmin(), crop_smooth.shape)
    return float(y0 + iy), float(x0 + ix)


# ── Capsid wall refinement ───────────────────────────────────────────────────

def _per_sector_radial_profile(
    img: np.ndarray,
    cy: float, cx: float,
    r_start: int, r_end: int,
    n_sectors: int = 8,
    n_angles_per_sector: int = 12,
) -> np.ndarray:
    """Radial profile per angular sector. Returns (n_sectors, n_radii)."""
    h, w = img.shape
    profile = np.zeros((n_sectors, r_end - r_start))
    sector_width = 2 * np.pi / n_sectors
    for s in range(n_sectors):
        a0 = s * sector_width
        a1 = (s + 1) * sector_width
        angles = np.linspace(a0, a1, n_angles_per_sector, endpoint=False)
        for i, r in enumerate(range(r_start, r_end)):
            ys = np.clip((cy + r * np.sin(angles)).astype(int), 0, h - 1)
            xs = np.clip((cx + r * np.cos(angles)).astype(int), 0, w - 1)
            profile[s, i] = img[ys, xs].mean()
    return profile


def assess_capsid_quality(
    img_norm: np.ndarray,
    cy: float, cx: float,
    wall_r_px: float,
    expected_r_px: float,
) -> dict:
    """
    Compute quality metrics for a candidate capsid at refined wall radius.

    Returns dict with:
      contrast        : (ring_mean - exterior_mean), should be > 0 (bright ring on dark stain)
      interior_uniformity : 1 - (ring_std / ring_mean), higher = more uniform protein ring
      wall_radius_cv  : std/mean of per-sector wall radii — lower = more circular
    """
    h, w = img_norm.shape
    # wall_r_px now points to the dark stain-ring minimum (outer boundary).
    # The bright protein ring sits INSIDE that radius; the post-stain
    # background sits AT the wall radius and just past it. So:
    #   ring annulus    : inside the wall, where bright protein lives
    #   exterior annulus: at the wall, where the dark stain pool is
    # contrast = bright_protein - dark_stain stays positive for real BMVs.
    yy, xx = np.ogrid[:h, :w]
    rr = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)

    ring     = (rr >= 0.55 * wall_r_px) & (rr <= 0.85 * wall_r_px)
    exterior = (rr >= 0.95 * wall_r_px) & (rr <= 1.10 * wall_r_px)

    if ring.sum() < 10 or exterior.sum() < 10:
        return dict(contrast=np.nan, interior_uniformity=np.nan,
                    wall_radius_cv=np.nan, exterior_mean=np.nan)

    ring_vals = img_norm[ring]
    ext_vals = img_norm[exterior]
    ring_mean = float(ring_vals.mean())
    ring_std = float(ring_vals.std())
    ext_mean = float(ext_vals.mean())

    contrast = ring_mean - ext_mean
    uniformity = 1.0 - (ring_std / ring_mean) if ring_mean > 1e-6 else 0.0
    # Exterior must actually be dark (stain pool), not just darker-than-ring.
    # Stored so caller can filter on absolute exterior intensity.
    exterior_mean = ext_mean

    # Per-sector wall radius
    r_start = max(1, int(expected_r_px * 0.5))
    r_end = min(int(expected_r_px * 1.6) + 1,
                min(h, w) // 2,
                int(min(cy, h - cy, cx, w - cx)) - 1)
    if r_end - r_start < 4:
        return dict(contrast=contrast, interior_uniformity=uniformity,
                    wall_radius_cv=np.nan, exterior_mean=exterior_mean)

    sector_prof = _per_sector_radial_profile(img_norm, cy, cx, r_start, r_end)
    sector_smooth = gaussian_filter1d(sector_prof, sigma=2.0, axis=1)
    sector_grad = np.gradient(sector_smooth, axis=1)
    sector_walls = r_start + sector_grad.argmin(axis=1)
    cv = float(sector_walls.std() / max(sector_walls.mean(), 1e-6))
    return dict(contrast=contrast, interior_uniformity=uniformity,
                wall_radius_cv=cv, exterior_mean=exterior_mean)


def refine_capsid_wall(
    img_norm: np.ndarray,
    cy: float,
    cx: float,
    expected_r_px: float,
    nm_per_px: float,
    smooth_sigma: float | None = None,
    r_start_frac: float = PROFILE_R_START_FRAC,
    r_end_frac: float = PROFILE_R_END_FRAC,
) -> float | None:
    """
    Find capsid outer wall as the inflection point of the radial profile.

    For uranyl-stained BMV capsids the radial profile is bright interior →
    monotonic descent → dark stain exterior (no clear minimum-then-recovery
    like the VLP/gold case). The wall is therefore the steepest descending
    point on the slope.

    Search is constrained to [r_start_frac, r_end_frac] × expected_r_px so
    that internal texture (RNA core, stain pools) and far-field background
    cannot win.
    """
    h, w = img_norm.shape
    r_start = max(1, int(expected_r_px * r_start_frac))
    r_end = int(expected_r_px * r_end_frac) + 1
    r_end = min(r_end, min(h, w) // 2,
                int(min(cy, h - cy, cx, w - cx)) - 1)
    if r_end - r_start < 4:
        return None

    profile = radial_profile(img_norm, cy, cx, r_start, r_end)
    sigma_px = smooth_sigma if smooth_sigma is not None else max(0.5, PROFILE_SMOOTH_NM / nm_per_px)
    smoothed = gaussian_filter1d(profile, sigma=sigma_px)

    # BMV donut profile: dark interior → bright protein-ring peak → dark
    # stain-pool minimum → background. We pick the wall as the radius
    # where intensity has fallen to a fraction WALL_DESCENT_FRAC of the way
    # from the peak down to the post-peak minimum. 0.5 = half-max
    # (matches manual outer-protein-ring measurement). Larger values move
    # the wall outward toward the dark-stain trough.
    peak_i = int(np.argmax(smoothed))
    if peak_i >= len(smoothed) - 2:
        return None

    descent = smoothed[peak_i:]
    post_min_i = int(np.argmin(descent))
    if post_min_i <= 0:
        return None
    peak_v = float(smoothed[peak_i])
    min_v = float(descent[post_min_i])
    if peak_v - min_v < 0.02:
        return None

    target = peak_v - WALL_DESCENT_FRAC * (peak_v - min_v)
    edge_local = None
    for j in range(1, post_min_i + 1):
        if descent[j] <= target:
            v0, v1 = descent[j - 1], descent[j]
            edge_local = (j - 1) + (v0 - target) / (v0 - v1) if v0 != v1 else j
            break
    if edge_local is None:
        return None
    return float(r_start + peak_i + edge_local)


# ── Output ───────────────────────────────────────────────────────────────────

def save_overlay(
    img_norm: np.ndarray,
    df: pd.DataFrame,
    nm_per_px: float,
    out_path: Path,
    title: str = "",
    highlight_indices: set | None = None,
    show_flagged: bool = False,
) -> None:
    """Cyan circle = refined capsid wall; magenta = blob-detected scale (fallback)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    axes[0].imshow(img_norm, cmap="gray", interpolation="none")
    axes[0].set_title("Image")
    axes[0].axis("off")

    axes[1].imshow(img_norm, cmap="gray", interpolation="none")
    reliable_col = df["is_reliable"] if "is_reliable" in df.columns else pd.Series(True, index=df.index)
    for idx, row in df.iterrows():
        cx, cy = row["centroid_x"], row["centroid_y"]
        highlighted = highlight_indices is not None and idx in highlight_indices
        is_reliable = bool(reliable_col.loc[idx])

        if not np.isnan(row["capsid_diameter_nm"]):
            if not (is_reliable or highlighted or show_flagged):
                continue
            r_px = row["capsid_diameter_nm"] / 2 / nm_per_px
            if highlighted:
                color, lw, alpha = "orange", 1.0, 1.0
            elif is_reliable:
                color, lw, alpha = "cyan", 0.6, 1.0
            else:
                color, lw, alpha = "tomato", 0.5, 0.85
        else:
            if not show_flagged:
                continue
            r_px = row["blob_radius_px"]
            color, lw, alpha = "magenta", 0.5, 0.85

        axes[1].add_patch(mpatches.Circle(
            (cx, cy), r_px,
            linewidth=lw, edgecolor=color, facecolor="none", alpha=alpha
        ))
        axes[1].plot([cx, cx + r_px], [cy, cy], color=color, linewidth=lw * 0.6, alpha=alpha)

        # Only label reliable + highlighted to keep overlay readable
        if is_reliable or highlighted:
            axes[1].text(
                cx, cy - r_px - 2, str(idx),
                color=color,
                fontsize=5 if highlighted else 3.5,
                fontweight="bold" if highlighted else "normal",
                ha="center", va="bottom",
                bbox=dict(boxstyle="round,pad=0.1", facecolor="black", alpha=0.6, linewidth=0),
            )

    n_total = len(df)
    n_reliable = int(reliable_col.sum())
    axes[1].set_title(f"n={n_total} candidates  ({n_reliable} reliable)")
    axes[1].axis("off")
    axes[1].legend(handles=[
        mpatches.Patch(color="cyan", label="Reliable capsid"),
        mpatches.Patch(color="tomato", label="Flagged (failed quality)"),
        mpatches.Patch(color="magenta", label="No wall fit"),
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
    expected_diameter_nm: float,
    search_factor: float,
    out_path: Path,
    indices: list[int] | None = None,
) -> None:
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

    expected_r_px = expected_diameter_nm / 2 / nm_per_px
    r_start = 1  # debug plot shows the full profile from center for context
    r_end = int(expected_r_px * search_factor) + 1

    for ax, (idx, row) in zip(axes, sample.iterrows()):
        cy, cx = row["centroid_y"], row["centroid_x"]
        profile = radial_profile(img_norm, cy, cx, r_start, r_end)
        smoothed = gaussian_filter1d(profile, sigma=3.0)
        radii_nm = (r_start + np.arange(len(profile))) * nm_per_px

        ax.plot(radii_nm, profile, color="lightsteelblue", linewidth=1, label="raw")
        ax.plot(radii_nm, smoothed, color="steelblue", linewidth=1.5, label="smoothed")
        ax.axvline(expected_diameter_nm / 2, color="gray", linewidth=1, linestyle=":",
                   label=f"expected r = {expected_diameter_nm/2:.1f} nm")
        if not np.isnan(row["capsid_diameter_nm"]):
            ax.axvline(row["capsid_diameter_nm"] / 2, color="cyan", linewidth=1.5,
                       label=f"refined r = {row['capsid_diameter_nm']/2:.1f} nm")
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


def save_histogram(df_all: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    data = df_all["capsid_diameter_nm"].dropna()
    if len(data) == 0:
        ax.set_title("Capsid diameter (no data)")
    else:
        bins = np.arange(np.floor(data.min()), np.ceil(data.max()) + 1, 1)
        ax.hist(data, bins=bins, edgecolor="white", linewidth=0.5, color="steelblue")
        ax.axvline(data.mean(), color="tomato", linewidth=2,
                   label=f"mean = {data.mean():.1f} ± {data.std():.1f} nm  (n={len(data)})")
        ax.set_xlabel("Capsid diameter (nm)")
        ax.set_ylabel("Count")
        ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Histogram  → {out_path}")


# ── Standalone progress GUI (daemon thread, closeable) ──────────────────────

class ProgressGUI:
    """
    Tkinter window that shows current progress in a separate, non-blocking
    daemon thread. User can close the window without killing the program;
    subsequent updates from the worker just become no-ops.

    On macOS, Tk normally insists on running in the main thread, but a
    small dedicated daemon thread works for a lightweight read-only window.
    If tkinter import or Tk() fails for any reason, the GUI silently
    becomes a no-op so the rest of the run is unaffected.
    """

    def __init__(self, total: int, title: str = "BMV measurement"):
        self.total = total
        self.q: queue.Queue = queue.Queue()
        self._closed = threading.Event()
        self._ready = threading.Event()
        self._thread = threading.Thread(
            target=self._run, args=(title,), daemon=True)
        self._thread.start()
        self._ready.wait(timeout=2.0)

    def _run(self, title: str) -> None:
        try:
            import tkinter as tk
            from tkinter import ttk
        except Exception:
            self._ready.set()
            return
        try:
            root = tk.Tk()
            root.title(title)
            root.geometry("420x140")
            self._root = root

            self._label = tk.Label(root, text="starting…",
                                   font=("Helvetica", 12))
            self._label.pack(pady=(18, 6))
            self._sub = tk.Label(root, text=f"0 / {self.total}",
                                 font=("Helvetica", 10), fg="#666")
            self._sub.pack()
            self._bar = ttk.Progressbar(root, length=380, mode="determinate",
                                        maximum=self.total)
            self._bar.pack(pady=10)
            tk.Button(root, text="Close", command=self._close).pack()
            root.protocol("WM_DELETE_WINDOW", self._close)

            def poll():
                try:
                    while True:
                        msg = self.q.get_nowait()
                        if isinstance(msg, tuple):
                            done, current = msg
                            self._bar["value"] = done
                            self._sub.config(text=f"{done} / {self.total}")
                            self._label.config(text=current)
                        else:
                            self._label.config(text=str(msg))
                except queue.Empty:
                    pass
                if not self._closed.is_set():
                    root.after(150, poll)

            self._ready.set()
            root.after(150, poll)
            root.mainloop()
        except Exception:
            self._ready.set()

    def _close(self) -> None:
        self._closed.set()
        try:
            self._root.destroy()
        except Exception:
            pass

    def update(self, done: int, current: str) -> None:
        if self._closed.is_set():
            return
        try:
            self.q.put_nowait((done, current))
        except Exception:
            pass


# ── Per-image worker ─────────────────────────────────────────────────────────

def process_image(
    path: Path,
    out_dir: Path,
    expected_nm: float,
    save_debug: bool,
    debug_indices: list[int] | None,
    show_flagged: bool,
) -> tuple[pd.DataFrame, str]:
    """
    Run the full pipeline on a single image. Top-level so it can be pickled
    by ProcessPoolExecutor. Returns (per-image dataframe, summary log lines).
    """
    log = [f"\n{path.name}"]
    img, nm_per_px = load_dm(path)
    log.append(f"  {img.shape[0]}×{img.shape[1]} px  |  {nm_per_px:.4f} nm/px")
    expected_r_px = expected_nm / 2 / nm_per_px

    img_norm = normalise(img)
    # Detection-only image: smoothed at the capsomere scale so high-mag
    # images don't have sub-capsid texture confusing Hough / center-snap.
    detect_sigma_px = max(0.5, DETECTION_PRESMOOTH_NM / nm_per_px)
    img_for_detect = gaussian_filter(img_norm, sigma=detect_sigma_px)

    df = detect_capsid_centers(img_for_detect, nm_per_px, expected_nm)
    log.append(f"  {len(df)} candidate centers detected (Hough, presmooth {detect_sigma_px:.1f}px)")

    # Center refinement (snap to local intensity min if it lowers wall_cv)
    # Center snap uses pre-smoothed image; wall fit uses raw img_norm.
    refined_y, refined_x = [], []
    for _, row in df.iterrows():
        cy0, cx0 = row["centroid_y"], row["centroid_x"]
        cy1, cx1 = refine_center(img_for_detect, cy0, cx0, expected_r_px)
        r0 = refine_capsid_wall(img_norm, cy0, cx0, expected_r_px, nm_per_px)
        r1 = refine_capsid_wall(img_norm, cy1, cx1, expected_r_px, nm_per_px)
        if r0 is not None and r1 is not None:
            cv0 = assess_capsid_quality(img_norm, cy0, cx0, r0, expected_r_px)["wall_radius_cv"]
            cv1 = assess_capsid_quality(img_norm, cy1, cx1, r1, expected_r_px)["wall_radius_cv"]
            if not np.isnan(cv1) and (np.isnan(cv0) or cv1 < cv0):
                refined_y.append(cy1); refined_x.append(cx1); continue
        refined_y.append(cy0); refined_x.append(cx0)
    df["centroid_y"] = refined_y
    df["centroid_x"] = refined_x

    # Wall fit + quality scoring
    capsid_d, contrasts, unifs, wall_cvs, ext_means = [], [], [], [], []
    for _, row in df.iterrows():
        r_px = refine_capsid_wall(img_norm, row["centroid_y"], row["centroid_x"],
                                  expected_r_px, nm_per_px)
        capsid_d.append(2 * r_px * nm_per_px if r_px else np.nan)
        if r_px is not None:
            q = assess_capsid_quality(img_norm, row["centroid_y"], row["centroid_x"],
                                      r_px, expected_r_px)
            contrasts.append(q["contrast"]); unifs.append(q["interior_uniformity"])
            wall_cvs.append(q["wall_radius_cv"]); ext_means.append(q["exterior_mean"])
        else:
            contrasts.append(np.nan); unifs.append(np.nan)
            wall_cvs.append(np.nan); ext_means.append(np.nan)
    df["capsid_diameter_nm"] = capsid_d
    df["contrast"] = contrasts
    df["interior_uniformity"] = unifs
    df["wall_radius_cv"] = wall_cvs
    df["exterior_mean"] = ext_means
    df = df.reset_index(drop=True)

    # Quality filtering
    lo = expected_nm * (1 - DIAM_TOL_FRAC)
    hi = expected_nm * (1 + DIAM_TOL_FRAC)
    df["is_reliable"] = (
        df["capsid_diameter_nm"].between(lo, hi)
        & (df["contrast"] >= MIN_CONTRAST)
        & (df["interior_uniformity"] >= MIN_UNIFORMITY)
        & (df["wall_radius_cv"] <= MAX_WALL_CV)
        & (df["exterior_mean"] <= MAX_EXTERIOR_MEAN)
    )

    # Overlap exclusion
    df["overlapping"] = False
    rel_idx = df.index[df["is_reliable"]].tolist()
    if len(rel_idx) >= 2:
        radii_px = (df.loc[rel_idx, "capsid_diameter_nm"] / 2 / nm_per_px).values
        pos = df.loc[rel_idx, ["centroid_x", "centroid_y"]].values
        tree = cKDTree(pos)
        for i, j in tree.query_pairs(2 * float(radii_px.max())):
            d = np.hypot(pos[i, 0] - pos[j, 0], pos[i, 1] - pos[j, 1])
            if d < (radii_px[i] + radii_px[j]) * OVERLAP_TOL:
                df.at[rel_idx[i], "overlapping"] = True
                df.at[rel_idx[j], "overlapping"] = True
        df.loc[df["overlapping"], "is_reliable"] = False

    n_refined = int(df["capsid_diameter_nm"].notna().sum())
    n_reliable = int(df["is_reliable"].sum())
    n_overlap = int(df["overlapping"].sum())
    log.append(f"  walls refined {n_refined}/{len(df)}, "
               f"{n_reliable} reliable, {n_overlap} overlap-dropped")
    if n_reliable:
        rel = df[df["is_reliable"]]["capsid_diameter_nm"]
        log.append(f"  mean {rel.mean():.1f} ± {rel.std():.1f} nm  (n={n_reliable})")

    df.insert(0, "file", path.name)

    # Overlay
    overlay_dir = out_dir / "overlays"
    overlay_dir.mkdir(parents=True, exist_ok=True)
    sample_idx = []
    if save_debug:
        if debug_indices:
            sample_idx = [i for i in debug_indices if i in df.index]
        else:
            sample_idx = (df[df["is_reliable"]]
                          .dropna(subset=["capsid_diameter_nm"])
                          .head(6).index.tolist())
    save_overlay(img_norm, df, nm_per_px, overlay_dir / (path.stem + "_overlay.png"),
                 title=path.name, highlight_indices=set(sample_idx) if sample_idx else None,
                 show_flagged=show_flagged)
    if save_debug:
        save_debug_profiles(img_norm, df, nm_per_px, expected_nm, PROFILE_R_END_FRAC,
                            overlay_dir / (path.stem + "_profiles.png"),
                            indices=sample_idx)

    return df, "\n".join(log)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("image_path", type=Path,
                        help="Image file (.dm3/.dm4) or folder of images")
    parser.add_argument("--pattern", default="*",
                        help="Glob when image_path is a folder")
    parser.add_argument("--out", type=Path, default=Path("results/bmv"),
                        help="Output directory (default: results/bmv/)")

    parser.add_argument("--expected-nm", type=float, default=28.0,
                        help="Expected capsid diameter in nm (default 28). "
                             "All other tuning lives as constants at the top of bmv_measure.py.")

    parser.add_argument("--debug", action="store_true",
                        help="Save radial profile plots for a sample of particles")
    parser.add_argument("--debug-indices", type=int, nargs="+", default=None)
    parser.add_argument("--show-flagged", action="store_true",
                        help="Show flagged (tomato) and no-fit (magenta) circles in the overlay. "
                             "By default only reliable cyan circles are drawn for cleaner inspection.")
    parser.add_argument("--workers", type=int, default=4,
                        help="Parallel worker processes for multi-image runs (default 4). "
                             "Set to 1 for serial / single-file runs.")
    parser.add_argument("--gui", action="store_true",
                        help="Show a separate Tk progress window. Closing the window does "
                             "not stop the run; the rest of the pipeline continues.")

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
    (args.out / "overlays").mkdir(exist_ok=True)

    all_rows = []

    gui = ProgressGUI(total=len(files)) if args.gui else None
    if gui:
        gui.update(0, "starting…")

    done = 0
    if args.workers > 1 and len(files) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {
                ex.submit(process_image, p, args.out, args.expected_nm,
                          args.debug, args.debug_indices, args.show_flagged): p
                for p in files
            }
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc="images", unit="img", ncols=80):
                df_i, log = fut.result()
                tqdm.write(log)
                all_rows.append(df_i)
                done += 1
                if gui:
                    gui.update(done, f"finished {futures[fut].name}")
    else:
        for path in tqdm(files, desc="images", unit="img", ncols=80):
            if gui:
                gui.update(done, f"processing {path.name}")
            df_i, log = process_image(path, args.out, args.expected_nm,
                                      args.debug, args.debug_indices,
                                      args.show_flagged)
            tqdm.write(log)
            all_rows.append(df_i)
            done += 1
            if gui:
                gui.update(done, f"finished {path.name}")

    if gui:
        gui.update(len(files), "done — close window")

    results = pd.concat(all_rows, ignore_index=True)

    csv_path = args.out / "bmv_measurements.csv"
    results.to_csv(csv_path, index=False)
    print(f"\n  CSV        → {csv_path}")

    hist_path = args.out / "bmv_histogram.png"
    save_histogram(results[results["is_reliable"]], hist_path)

    print("\n── Summary ──────────────────────────────────────────────")
    n_flagged = (~results["is_reliable"] & results["capsid_diameter_nm"].notna()).sum()
    reliable = results[results["is_reliable"]]["capsid_diameter_nm"].dropna()
    print(f"\nCapsid — reliable (n={len(reliable)}, {n_flagged} out-of-range flagged):")
    if len(reliable):
        print(f"  mean    {reliable.mean():.1f} nm")
        print(f"  std     {reliable.std():.1f} nm")
        print(f"  median  {reliable.median():.1f} nm")
        print(f"  range   {reliable.min():.1f} – {reliable.max():.1f} nm")

    # Auto-eval against benchmarks/bmv/ (if available). Wrapped in try/except
    # so this script remains usable as a single-file standalone — drop alone
    # next to images and it still measures.
    # Build a run-result dict so we can reuse the VLP summary writer + eval pipeline.
    try:
        from analysis.vlp_measure_v2 import (  # noqa: WPS433
            _per_image_summary, _overall_summary, _write_summary_md,
        )
        run_result = {
            "sample_type":    "BMV",
            "sample_name":    args.out.name,
            "script_version": "bmv_measure@1.0",
            "summary":        _overall_summary(results),
            "per_image":      _per_image_summary(results),
            "outputs":        {"run_dir":         str(args.out),
                               "csv_path":        str(csv_path),
                               "overlays_dir":    str(args.out / "overlays"),
                               "histograms_path": str(hist_path),
                               "summary_md":      ""},
        }
    except ImportError:
        run_result = None  # standalone single-file use; skip summary + eval

    # Auto-eval against benchmarks/bmv/ (if available).
    if run_result is not None:
        try:
            from analysis.eval import evaluate  # noqa: WPS433
            ev = evaluate(run_result, sample_type="BMV", write_report=True)
            run_result["eval"] = {
                "n_warns_total":  ev["n_warns_total"],
                "n_ref_runs":     ev["n_ref_runs"],
                "hand_vs_script": ev["hand_vs_script"],
                "report_path":    ev.get("report_path"),
            }
            print(f"\n  Eval → {ev.get('report_path', '(no report)')}")
            print(f"        warnings: {ev['n_warns_total']} across {ev['n_images']} image(s)")
            if ev["hand_vs_script"]:
                H = ev["hand_vs_script"]
                print(f"        hand vs script: Δ capsid {H['delta_capsid_nm']:+.2f} nm "
                      f"(hand n={H['n_hand_particles']})")
        except FileNotFoundError:
            run_result["eval"] = {"skipped": "no reference_runs.csv found"}
        except Exception as exc:  # pragma: no cover (defensive)
            run_result["eval"] = {"skipped": f"eval error: {exc}"}
            print(f"  (eval skipped: {exc})", file=sys.stderr)

        # SUMMARY.md last so the eval block can be folded in.
        summary_md_path = _write_summary_md(args.out, run_result)
        print(f"  SUMMARY.md → {summary_md_path}")


if __name__ == "__main__":
    main()
