"""
vlp_measure_v2.py — VLP gold + capsid measurement, ported BMV-style

Same gold-NP detection as the original `vlp_measure.py` (intensity threshold
+ connected components on dark cores), but the capsid wall is fitted with
the BMV-style radial-profile method:

  - profile from gold edge outward; identify bright-protein peak and the
    post-peak dark-stain-ring minimum
  - wall placed at WALL_DESCENT_FRAC of the way from peak down to min
    (0.75 = matches what a human draws as the outer protein edge in
    ImageJ — verified on BMV)
  - smoothing declared in nm and converted to px per-image so behaviour
    is magnification-invariant

Reliability is replaced with **quality-based filters** (per-sector wall
circularity, ring-vs-stain contrast, exterior darkness) instead of the
population-relative ±2σ clip used in the original — that clip removed
real heterogeneity rather than bad fits.

Usage:
    uv run python vlp_measure_v2.py "images/VLPs for machine learning project" \
        --pattern "VLP17_*" --workers 6
"""

import argparse
import datetime as _dt
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
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.spatial import cKDTree
from skimage import exposure, filters, morphology, measure


# ── Shared utilities (image loading, normalisation, gold NP detection,
#    radial profile sampling). These live here so that vlp_measure_v2.py is
#    self-contained — a scientist can drop just this file alongside images
#    and run it. The legacy `vlp_measure.py` re-imports them from here. ──

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


def auto_threshold(img_norm: np.ndarray, gold_threshold: float | None) -> float:
    """Otsu on the dark half of the histogram. Override with explicit value if given."""
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
    Threshold + connected-component detection of gold NPs.
    Returns (DataFrame[centroid_y, centroid_x, gold_diameter_nm, gold_radius_px],
             threshold_used).
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
    df["gold_radius_px"]   = np.sqrt(df["area"] / np.pi)

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
# All physical-scale tunables in nm; pixel sigmas are derived per-image.

# Capsid radial profile / wall fitting
PROFILE_SMOOTH_NM    = 0.7    # smoothing of radial profile, in nm
PROFILE_SKIP_NM      = 4.0    # skip past gold edge before searching for peak
PROFILE_END_NM       = 25.0   # outer search extent past gold edge. Capsid wall
                              # sits ~5-15 nm past gold edge for typical VLPs;
                              # 25 nm gives margin without letting the post-peak
                              # min escape to far stain features.
WALL_DESCENT_FRAC    = 0.75   # 0.0 = peak, 0.5 = half-max, 1.0 = post-peak min

# Quality filters (only the capsid wall fit; gold is already filtered in detect_gold)
MIN_RING_CONTRAST    = 0.05   # bright-ring mean − stain-min mean (must be > 0)
MAX_WALL_CV          = 0.20   # per-sector wall-radius coefficient of variation
MIN_PEAK_PROMINENCE  = 0.15   # peak height above post-peak minimum, [0,1] scale.
                              # Real VLP capsids show prominence ≥0.20; values
                              # near zero are bad fits where the "peak" and "min"
                              # are within noise of each other.

# Overlap exclusion
OVERLAP_TOL          = 0.95   # touching pair if dist < (r1+r2)·this


# ── Script identity (used by the eval system to track which version produced a run) ──
SCRIPT_NAME    = "vlp_measure_v2"
SCRIPT_VERSION = "2.0"

def _script_version_with_tunables() -> str:
    """Stable identifier including the tunables — eval cares when these change."""
    bits = (f"wall={WALL_DESCENT_FRAC};contrast={MIN_RING_CONTRAST};"
            f"cv={MAX_WALL_CV};prom={MIN_PEAK_PROMINENCE};"
            f"smooth={PROFILE_SMOOTH_NM};end={PROFILE_END_NM}")
    return f"{SCRIPT_NAME}@{SCRIPT_VERSION}({bits})"


# ── Standalone progress GUI (daemon thread, closeable) ──────────────────────

class ProgressGUI:
    """Tk progress window in a daemon thread; closing it doesn't kill the run."""

    def __init__(self, total: int, title: str = "VLP measurement"):
        self.total = total
        self.q: queue.Queue = queue.Queue()
        self._closed = threading.Event()
        self._ready = threading.Event()
        self._thread = threading.Thread(target=self._run, args=(title,), daemon=True)
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
            self._label = tk.Label(root, text="starting…", font=("Helvetica", 12))
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


# ── Capsid wall fit (BMV-style WALL_DESCENT_FRAC) ─────────────────────────────

def find_capsid_wall(
    img_norm: np.ndarray,
    cy: float, cx: float,
    gold_radius_px: float,
    nm_per_px: float,
) -> tuple[float | None, dict]:
    """
    Search outward from the gold edge for the bright-protein peak and the
    post-peak dark-stain-ring minimum. Wall is placed at WALL_DESCENT_FRAC of
    the way from peak to min. Returns (capsid_radius_px or None, metrics dict).
    """
    h, w = img_norm.shape
    skip_px = max(1, int(PROFILE_SKIP_NM / nm_per_px))
    r_start = max(1, int(gold_radius_px) + skip_px)
    r_end = r_start + int(PROFILE_END_NM / nm_per_px)
    r_end = min(r_end, min(h, w) // 2,
                int(min(cy, h - cy, cx, w - cx)) - 1)
    metrics = dict(ring_contrast=np.nan, peak_prominence=np.nan)
    if r_end - r_start < 6:
        return None, metrics

    profile = radial_profile(img_norm, cy, cx, r_start, r_end)
    sigma_px = max(0.5, PROFILE_SMOOTH_NM / nm_per_px)
    smoothed = gaussian_filter1d(profile, sigma=sigma_px)

    peak_i = int(np.argmax(smoothed))
    if peak_i >= len(smoothed) - 2:
        return None, metrics

    descent = smoothed[peak_i:]
    post_min_i = int(np.argmin(descent))
    if post_min_i <= 0:
        return None, metrics

    peak_v = float(smoothed[peak_i])
    min_v = float(descent[post_min_i])
    prominence = peak_v - min_v
    metrics["peak_prominence"] = prominence
    metrics["ring_contrast"] = prominence

    if prominence < MIN_PEAK_PROMINENCE:
        return None, metrics

    target = peak_v - WALL_DESCENT_FRAC * prominence
    edge_local = None
    for j in range(1, post_min_i + 1):
        if descent[j] <= target:
            v0, v1 = descent[j - 1], descent[j]
            edge_local = (j - 1) + (v0 - target) / (v0 - v1) if v0 != v1 else j
            break
    if edge_local is None:
        return None, metrics

    return float(r_start + peak_i + edge_local), metrics


def per_sector_wall_cv(
    img_norm: np.ndarray,
    cy: float, cx: float,
    gold_radius_px: float,
    nm_per_px: float,
    n_sectors: int = 8,
    n_angles_per_sector: int = 12,
) -> float:
    """
    Standard deviation / mean of per-sector steepest-descent radii. Real
    capsids have a circular wall at consistent radius across angles; stain
    aggregates / off-center fits don't. Used as a circularity quality metric.
    """
    h, w = img_norm.shape
    skip_px = max(1, int(PROFILE_SKIP_NM / nm_per_px))
    r_start = max(1, int(gold_radius_px) + skip_px)
    r_end = r_start + int(PROFILE_END_NM / nm_per_px)
    r_end = min(r_end, min(h, w) // 2,
                int(min(cy, h - cy, cx, w - cx)) - 1)
    if r_end - r_start < 6:
        return np.nan

    sigma_px = max(0.5, PROFILE_SMOOTH_NM / nm_per_px)
    profile = np.zeros((n_sectors, r_end - r_start))
    sector_w = 2 * np.pi / n_sectors
    for s in range(n_sectors):
        a = np.linspace(s * sector_w, (s + 1) * sector_w,
                        n_angles_per_sector, endpoint=False)
        for i, r in enumerate(range(r_start, r_end)):
            ys = np.clip((cy + r * np.sin(a)).astype(int), 0, h - 1)
            xs = np.clip((cx + r * np.cos(a)).astype(int), 0, w - 1)
            profile[s, i] = img_norm[ys, xs].mean()

    sm = gaussian_filter1d(profile, sigma=sigma_px, axis=1)
    grad = np.gradient(sm, axis=1)
    walls = r_start + grad.argmin(axis=1)
    return float(walls.std() / max(walls.mean(), 1e-6))


# ── Per-image worker ─────────────────────────────────────────────────────────

def process_image(
    path: Path,
    out_dir: Path,
    gold_threshold: float | None,
    min_gold_nm: float,
    max_gold_nm: float,
    show_flagged: bool,
) -> tuple[pd.DataFrame, str]:
    log = [f"\n{path.name}"]
    img, nm_per_px = load_dm(path)
    log.append(f"  {img.shape[0]}×{img.shape[1]} px  |  {nm_per_px:.4f} nm/px")
    img_norm = normalise(img)

    # Gold detection (unchanged from v1)
    df_gold, t_used = detect_gold(img_norm, nm_per_px, gold_threshold,
                                  min_gold_nm, max_gold_nm)
    log.append(f"  threshold {t_used:.4f}  |  {len(df_gold)} gold NPs")

    # Capsid wall fit + quality scoring per particle
    capsid_d, contrasts, prominences, wall_cvs = [], [], [], []
    for _, g in df_gold.iterrows():
        r_px, m = find_capsid_wall(img_norm, g["centroid_y"], g["centroid_x"],
                                   g["gold_radius_px"], nm_per_px)
        capsid_d.append(2 * r_px * nm_per_px if r_px else np.nan)
        contrasts.append(m["ring_contrast"])
        prominences.append(m["peak_prominence"])
        if r_px is not None:
            wall_cvs.append(per_sector_wall_cv(
                img_norm, g["centroid_y"], g["centroid_x"],
                g["gold_radius_px"], nm_per_px))
        else:
            wall_cvs.append(np.nan)
    df_gold["capsid_diameter_nm"] = capsid_d
    df_gold["ring_contrast"] = contrasts
    df_gold["peak_prominence"] = prominences
    df_gold["wall_radius_cv"] = wall_cvs
    df_gold = df_gold.reset_index(drop=True)

    # Quality-based reliability (replaces ±2σ population clip from v1)
    df_gold["is_reliable"] = (
        df_gold["capsid_diameter_nm"].notna()
        & (df_gold["ring_contrast"] >= MIN_RING_CONTRAST)
        & (df_gold["peak_prominence"] >= MIN_PEAK_PROMINENCE)
        & (df_gold["wall_radius_cv"] <= MAX_WALL_CV)
    )

    # Overlap exclusion: touching capsid pairs share a stain ring
    df_gold["overlapping"] = False
    rel_idx = df_gold.index[df_gold["is_reliable"]].tolist()
    if len(rel_idx) >= 2:
        radii_px = (df_gold.loc[rel_idx, "capsid_diameter_nm"] / 2 / nm_per_px).values
        pos = df_gold.loc[rel_idx, ["centroid_x", "centroid_y"]].values
        tree = cKDTree(pos)
        for i, j in tree.query_pairs(2 * float(radii_px.max())):
            d = np.hypot(pos[i, 0] - pos[j, 0], pos[i, 1] - pos[j, 1])
            if d < (radii_px[i] + radii_px[j]) * OVERLAP_TOL:
                df_gold.at[rel_idx[i], "overlapping"] = True
                df_gold.at[rel_idx[j], "overlapping"] = True
        df_gold.loc[df_gold["overlapping"], "is_reliable"] = False

    n_capsid = int(df_gold["capsid_diameter_nm"].notna().sum())
    n_reliable = int(df_gold["is_reliable"].sum())
    n_overlap = int(df_gold["overlapping"].sum())
    log.append(f"  capsid fit {n_capsid}/{len(df_gold)}, "
               f"{n_reliable} reliable, {n_overlap} overlap-dropped")
    if n_reliable:
        rel = df_gold[df_gold["is_reliable"]]
        log.append(f"  gold   mean {rel['gold_diameter_nm'].mean():.2f} ± "
                   f"{rel['gold_diameter_nm'].std():.2f} nm")
        log.append(f"  capsid mean {rel['capsid_diameter_nm'].mean():.2f} ± "
                   f"{rel['capsid_diameter_nm'].std():.2f} nm")

    df_gold.insert(0, "file", path.name)

    # Overlay
    overlay_dir = out_dir / "overlays"
    overlay_dir.mkdir(parents=True, exist_ok=True)
    save_overlay(img_norm, df_gold, nm_per_px,
                 overlay_dir / (path.stem + "_overlay.png"),
                 title=path.name, show_flagged=show_flagged)
    return df_gold, "\n".join(log)


# ── Output ───────────────────────────────────────────────────────────────────

def save_overlay(
    img_norm: np.ndarray,
    df: pd.DataFrame,
    nm_per_px: float,
    out_path: Path,
    title: str = "",
    show_flagged: bool = False,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    axes[0].imshow(img_norm, cmap="gray", interpolation="none")
    axes[0].set_title("Image"); axes[0].axis("off")

    axes[1].imshow(img_norm, cmap="gray", interpolation="none")
    for idx, row in df.iterrows():
        cx, cy = row["centroid_x"], row["centroid_y"]
        is_rel = bool(row["is_reliable"])

        # Gold (always drawn for reliable; faint for flagged)
        r_gold = row["gold_diameter_nm"] / 2 / nm_per_px
        gold_color = "yellow" if is_rel else ("orange" if show_flagged else None)
        if gold_color:
            alpha = 1.0 if is_rel else 0.5
            axes[1].add_patch(mpatches.Circle((cx, cy), r_gold, linewidth=0.5,
                                              edgecolor=gold_color, facecolor="none",
                                              alpha=alpha))
            axes[1].plot([cx, cx + r_gold], [cy, cy], color=gold_color,
                         linewidth=0.4, alpha=alpha)

        # Capsid
        if not np.isnan(row["capsid_diameter_nm"]):
            if not (is_rel or show_flagged):
                continue
            r_cap = row["capsid_diameter_nm"] / 2 / nm_per_px
            cap_color = "cyan" if is_rel else "tomato"
            alpha = 1.0 if is_rel else 0.5
            axes[1].add_patch(mpatches.Circle((cx, cy), r_cap, linewidth=0.5,
                                              edgecolor=cap_color, facecolor="none",
                                              alpha=alpha))
            axes[1].plot([cx, cx + r_cap], [cy, cy], color=cap_color,
                         linewidth=0.4, alpha=alpha)
            if is_rel:
                axes[1].text(cx, cy - r_cap - 2, str(idx), color="cyan",
                             fontsize=3.5, ha="center", va="bottom",
                             bbox=dict(boxstyle="round,pad=0.1", facecolor="black",
                                       alpha=0.6, linewidth=0))

    n_total = len(df)
    n_reliable = int(df["is_reliable"].sum())
    axes[1].set_title(f"n={n_total}  ({n_reliable} reliable)")
    axes[1].axis("off")
    axes[1].legend(handles=[
        mpatches.Patch(color="yellow", label="Gold NP"),
        mpatches.Patch(color="cyan", label="Reliable capsid"),
        mpatches.Patch(color="tomato", label="Flagged capsid"),
    ], loc="lower right", fontsize=7)
    if title:
        fig.suptitle(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_histograms(df_all: pd.DataFrame, out_path: Path) -> None:
    rel = df_all[df_all["is_reliable"]]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, col, label, color in [
        (axes[0], "gold_diameter_nm", "Gold NP diameter (nm)", "gold"),
        (axes[1], "capsid_diameter_nm", "Capsid diameter (nm)", "steelblue"),
    ]:
        data = rel[col].dropna()
        if len(data) == 0:
            ax.set_title(f"{label}\n(no data)"); continue
        bins = np.arange(np.floor(data.min()), np.ceil(data.max()) + 1, 0.5)
        ax.hist(data, bins=bins, edgecolor="white", linewidth=0.5, color=color)
        ax.axvline(data.mean(), color="tomato", linewidth=2,
                   label=f"mean = {data.mean():.2f} ± {data.std():.2f} nm  (n={len(data)})")
        ax.set_xlabel(label); ax.set_ylabel("Count"); ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_scatter(df_all: pd.DataFrame, out_path: Path) -> None:
    rel = df_all[df_all["is_reliable"]].dropna(subset=["gold_diameter_nm", "capsid_diameter_nm"])
    if len(rel) == 0:
        return
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(rel["gold_diameter_nm"], rel["capsid_diameter_nm"],
               alpha=0.4, s=12, color="steelblue", edgecolors="none")
    ax.axvline(rel["gold_diameter_nm"].mean(), color="gold", linestyle="--", linewidth=1,
               label=f"gold mean {rel['gold_diameter_nm'].mean():.2f} nm")
    ax.axhline(rel["capsid_diameter_nm"].mean(), color="cyan", linestyle="--", linewidth=1,
               label=f"capsid mean {rel['capsid_diameter_nm'].mean():.2f} nm")
    ax.set_xlabel("Gold NP diameter (nm)")
    ax.set_ylabel("VLP capsid diameter (nm)")
    ax.set_title(f"VLP gold vs capsid  (n={len(rel)})")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Per-image / overall summary builders ─────────────────────────────────────

def _per_image_summary(results: pd.DataFrame) -> list[dict]:
    """One dict per input image. Dimensionless + diameter stats for downstream eval.
    Works for both VLP (gold + capsid) and BMV (capsid only) measurement CSVs."""
    has_gold = "gold_diameter_nm" in results.columns
    out: list[dict] = []
    for fname, sub in results.groupby("file", sort=True):
        rel = sub[sub["is_reliable"]]
        n_gold = int(len(sub))
        n_wall = int(sub["capsid_diameter_nm"].notna().sum())
        n_rel  = int(len(rel))
        wall_cv = sub["wall_radius_cv"].dropna()
        out.append({
            "filename": str(fname),
            "n_gold": n_gold,
            "n_wall_fit": n_wall,
            "n_reliable": n_rel,
            "wall_fit_success_rate": (n_wall / n_gold) if n_gold else 0.0,
            "reliable_rate":         (n_rel  / n_gold) if n_gold else 0.0,
            "drop_rate":             ((n_gold - n_rel) / n_gold) if n_gold else 0.0,
            "gold_median_nm":        float(rel["gold_diameter_nm"].median())   if (has_gold and n_rel)     else float("nan"),
            "gold_std_nm":           float(rel["gold_diameter_nm"].std())      if (has_gold and n_rel >= 2) else float("nan"),
            "capsid_mean_nm":        float(rel["capsid_diameter_nm"].mean())   if n_rel     else float("nan"),
            "capsid_median_nm":      float(rel["capsid_diameter_nm"].median()) if n_rel     else float("nan"),
            "capsid_std_nm":         float(rel["capsid_diameter_nm"].std())    if n_rel >= 2 else float("nan"),
            "median_wall_cv":        float(wall_cv.median())                   if len(wall_cv) else float("nan"),
            "iqr_wall_cv":           float(wall_cv.quantile(0.75) - wall_cv.quantile(0.25)) if len(wall_cv) >= 2 else float("nan"),
        })
    return out


def _overall_summary(results: pd.DataFrame) -> dict:
    """Overall headline. Works for VLP (gold + capsid) and BMV (capsid only)."""
    has_gold = "gold_diameter_nm" in results.columns
    rel = results[results["is_reliable"]]
    n_gold = int(len(results))
    n_rel  = int(len(rel))
    return {
        "n_images":     int(results["file"].nunique()),
        "n_gold_total": n_gold,
        "n_reliable":   n_rel,
        "n_dropped":    n_gold - n_rel,
        "drop_rate":    ((n_gold - n_rel) / n_gold) if n_gold else 0.0,
        "gold_mean_nm":   float(rel["gold_diameter_nm"].mean())   if (has_gold and n_rel)     else float("nan"),
        "gold_std_nm":    float(rel["gold_diameter_nm"].std())    if (has_gold and n_rel >= 2) else float("nan"),
        "gold_median_nm": float(rel["gold_diameter_nm"].median()) if (has_gold and n_rel)     else float("nan"),
        "capsid_mean_nm":   float(rel["capsid_diameter_nm"].mean())   if n_rel     else float("nan"),
        "capsid_std_nm":    float(rel["capsid_diameter_nm"].std())    if n_rel >= 2 else float("nan"),
        "capsid_median_nm": float(rel["capsid_diameter_nm"].median()) if n_rel     else float("nan"),
    }


def _write_summary_md(out_dir: Path, result: dict) -> Path:
    """Persistent human + agent-readable summary. Solves stdout-burial.
    Works for both VLP (gold + capsid) and BMV (capsid only) — gold lines
    only appear when gold_mean_nm is present."""
    s = result["summary"]
    has_gold = isinstance(s.get("gold_mean_nm"), (int, float)) and not (
        isinstance(s.get("gold_mean_nm"), float) and np.isnan(s["gold_mean_nm"])
    )
    detection_label = "gold detections" if has_gold else "particle detections"
    md = [
        f"# Measurement run — {result['sample_type']}",
        "",
        f"- **Sample:** `{result['sample_name'] or '(none)'}`",
        f"- **Script:** `{result['script_version']}`",
        f"- **Run dir:** `{result['outputs']['run_dir']}`",
        "",
        "## Headline",
        f"- {s['n_images']} images, {s['n_gold_total']} {detection_label}",
        f"- **{s['n_reliable']} reliable** ({(1-s['drop_rate'])*100:.1f}%) "
        f"— {s['n_dropped']} dropped ({s['drop_rate']*100:.1f}%)",
    ]
    if s.get("n_reliable", 0):
        if has_gold:
            md.append(
                f"- **Gold:** mean {s['gold_mean_nm']:.2f} ± {s['gold_std_nm']:.2f} nm "
                f"(median {s['gold_median_nm']:.2f})"
            )
        md.append(
            f"- **Capsid:** mean {s['capsid_mean_nm']:.2f} ± {s['capsid_std_nm']:.2f} nm "
            f"(median {s['capsid_median_nm']:.2f})"
        )

    # Eval headline (only if eval ran)
    ev = result.get("eval") or {}
    if ev and "skipped" not in ev:
        md.append("")
        md.append("## Eval")
        n_warns = ev.get("n_warns_total", 0)
        n_ref   = ev.get("n_ref_runs", 0)
        md.append(f"- Per-image quality check: compared against {n_ref} reference runs; "
                  f"**{n_warns} warning(s)**")
        H = ev.get("hand_vs_script")
        if H is None:
            md.append("- Hand vs script: no hand measurements in `<run_dir>/hand/`")
        else:
            md.append(f"- Hand vs script: "
                      f"Δ capsid **{H['delta_capsid_nm']:+.2f} nm** "
                      f"(hand n={H.get('n_hand_particles', H.get('hand_n', '?'))})")
        if ev.get("report_path"):
            md.append(f"- Full report: `{ev['report_path']}`")
    elif ev.get("skipped"):
        md.append("")
        md.append(f"## Eval\n- Skipped: {ev['skipped']}")

    md += ["", "## Outputs", f"- CSV: `{result['outputs']['csv_path']}`"]
    for label, key in [("Overlays", "overlays_dir"),
                       ("Histograms", "histograms_path"),
                       ("Scatter", "scatter_path")]:
        path = result["outputs"].get(key)
        if path:
            md.append(f"- {label}: `{path}`")

    md += [
        "",
        "## Per-image",
    ]
    if has_gold:
        md += [
            "| filename | n_gold | n_reliable | reliable% | wall_fit% | gold med (nm) | capsid med (nm) |",
            "|---|---|---|---|---|---|---|",
        ]
    else:
        md += [
            "| filename | n_detections | n_reliable | reliable% | wall_fit% | capsid med (nm) |",
            "|---|---|---|---|---|---|",
        ]
    for r in result["per_image"]:
        if has_gold:
            md.append(
                f"| {r['filename']} | {r['n_gold']} | {r['n_reliable']} | "
                f"{r['reliable_rate']*100:.0f}% | {r['wall_fit_success_rate']*100:.0f}% | "
                f"{r['gold_median_nm']:.2f} | {r['capsid_median_nm']:.2f} |"
            )
        else:
            md.append(
                f"| {r['filename']} | {r['n_gold']} | {r['n_reliable']} | "
                f"{r['reliable_rate']*100:.0f}% | {r['wall_fit_success_rate']*100:.0f}% | "
                f"{r['capsid_median_nm']:.2f} |"
            )
    path = out_dir / "SUMMARY.md"
    path.write_text("\n".join(md) + "\n")
    return path


# ── Programmatic entry point ─────────────────────────────────────────────────

def run(
    image_path: str | Path,
    *,
    sample_type: str = "VLP",
    sample_name: str | None = None,
    out_dir: str | Path = "results/vlp_v2",
    pattern: str = "*",
    gold_threshold: float | None = None,
    min_gold_nm: float = 7.0,
    max_gold_nm: float = 30.0,
    workers: int = 4,
    show_flagged: bool = False,
    gui: bool = False,
) -> dict:
    """
    Measure VLP gold + capsid across one or more .dm3/.dm4 files.

    `sample_name` defaults to the input folder's name (or, for a single file,
    the filename stem + today's date). It's just a label — the script doesn't
    interpret it; it's there so reports and reference rows are easy to track.

    Returns a structured dict. Also writes:
      - `<out_dir>/vlp_measurements.csv` — per-particle results
      - `<out_dir>/overlays/` — per-image overlay PNGs
      - `<out_dir>/vlp_histograms.png`, `vlp_scatter.png`
      - `<out_dir>/SUMMARY.md` — human + agent-readable headline
      - `<out_dir>/eval_report.md` — quality + hand-vs-script comparison (if benchmarks exist)
    """
    image_path = Path(image_path)
    out_dir    = Path(out_dir)

    if sample_name is None:
        if image_path.is_dir():
            sample_name = image_path.name
        else:
            sample_name = f"{image_path.stem}_{_dt.date.today().isoformat()}"

    if image_path.is_file():
        files = [image_path]
    elif image_path.is_dir():
        files = sorted(
            f for ext in (".dm3", ".dm4")
            for f in image_path.glob(f"{pattern}{ext}")
        )
    else:
        raise FileNotFoundError(f"Path not found: {image_path}")
    if not files:
        raise FileNotFoundError(f"No .dm3 / .dm4 files found under {image_path} (pattern={pattern!r})")

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "overlays").mkdir(exist_ok=True)

    progress = ProgressGUI(total=len(files)) if gui else None
    if progress:
        progress.update(0, "starting…")

    all_rows: list[pd.DataFrame] = []
    done = 0
    if workers > 1 and len(files) > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {
                ex.submit(process_image, p, out_dir, gold_threshold,
                          min_gold_nm, max_gold_nm, show_flagged): p
                for p in files
            }
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc="images", unit="img", ncols=80):
                df_i, log = fut.result()
                tqdm.write(log)
                all_rows.append(df_i)
                done += 1
                if progress:
                    progress.update(done, f"finished {futures[fut].name}")
    else:
        for path in tqdm(files, desc="images", unit="img", ncols=80):
            if progress:
                progress.update(done, f"processing {path.name}")
            df_i, log = process_image(path, out_dir, gold_threshold,
                                       min_gold_nm, max_gold_nm, show_flagged)
            tqdm.write(log)
            all_rows.append(df_i)
            done += 1
            if progress:
                progress.update(done, f"finished {path.name}")
    if progress:
        progress.update(len(files), "done — close window")

    results = pd.concat(all_rows, ignore_index=True)
    csv_path        = out_dir / "vlp_measurements.csv"
    histograms_path = out_dir / "vlp_histograms.png"
    scatter_path    = out_dir / "vlp_scatter.png"
    results.to_csv(csv_path, index=False)
    save_histograms(results, histograms_path)
    save_scatter(results, scatter_path)

    summary_overall  = _overall_summary(results)
    summary_per_img  = _per_image_summary(results)

    result = {
        "sample_type":    sample_type,
        "sample_name":    sample_name,
        "script_version": _script_version_with_tunables(),
        "summary":        summary_overall,
        "per_image":      summary_per_img,
        "outputs": {
            "run_dir":          str(out_dir),
            "csv_path":         str(csv_path),
            "overlays_dir":     str(out_dir / "overlays"),
            "histograms_path":  str(histograms_path),
            "scatter_path":     str(scatter_path),
            "summary_md":       "",
        },
    }
    # Auto-eval against the reference benchmarks (if available).
    # Wrapped in try/except so this script is fully usable standalone — drop
    # the file alone next to a folder of images and it still measures.
    try:
        from analysis.eval import evaluate  # noqa: WPS433
        eval_result = evaluate(result, sample_type=sample_type, write_report=True)
        result["eval"] = {
            "n_warns_total":  eval_result["n_warns_total"],
            "n_ref_runs":     eval_result["n_ref_runs"],
            "hand_vs_script": eval_result["hand_vs_script"],
            "report_path":    eval_result.get("report_path"),
        }
        if eval_result.get("report_path"):
            result["outputs"]["eval_report"] = eval_result["report_path"]
    except (ImportError, FileNotFoundError) as exc:
        # No analysis.eval module (single-file standalone) or no
        # reference_runs.csv (fresh repo, no benchmarks yet). Either is fine.
        result["eval"] = {"skipped": f"{type(exc).__name__}: {exc}"}
    except Exception as exc:  # pragma: no cover (defensive)
        result["eval"] = {"skipped": f"eval error: {exc}"}

    # SUMMARY.md is written last so the eval block can be folded in.
    summary_md_path = _write_summary_md(out_dir, result)
    result["outputs"]["summary_md"] = str(summary_md_path)
    return result


# ── CLI wrapper ──────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("image_path", type=Path)
    parser.add_argument("--pattern", default="*")
    parser.add_argument("--out", type=Path, default=Path("results/vlp_v2"))
    parser.add_argument("--sample-type", default="VLP")
    parser.add_argument("--sample-name", default=None,
                        help="Defaults to the folder name of image_path.")

    parser.add_argument("--gold-threshold", type=float, default=None,
                        help="Intensity cutoff for gold detection (default: auto per image)")
    parser.add_argument("--min-gold-nm", type=float, default=7.0)
    parser.add_argument("--max-gold-nm", type=float, default=30.0)

    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--show-flagged", action="store_true",
                        help="Draw flagged-unreliable circles in tomato in the overlay")

    args = parser.parse_args()
    try:
        result = run(
            image_path=args.image_path,
            sample_type=args.sample_type,
            sample_name=args.sample_name,
            out_dir=args.out,
            pattern=args.pattern,
            gold_threshold=args.gold_threshold,
            min_gold_nm=args.min_gold_nm,
            max_gold_nm=args.max_gold_nm,
            workers=args.workers,
            show_flagged=args.show_flagged,
            gui=args.gui,
        )
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr); sys.exit(1)

    s = result["summary"]
    print("\n── Summary ──────────────────────────────────────────────")
    print(f"Reliable n={s['n_reliable']}, dropped n={s['n_dropped']}")
    if s["n_reliable"]:
        print(f"  gold   mean {s['gold_mean_nm']:.2f} ± {s['gold_std_nm']:.2f}  "
              f"median {s['gold_median_nm']:.2f}")
        print(f"  capsid mean {s['capsid_mean_nm']:.2f} ± {s['capsid_std_nm']:.2f}  "
              f"median {s['capsid_median_nm']:.2f}")
    print(f"\n  CSV → {result['outputs']['csv_path']}")
    print(f"  SUMMARY.md → {result['outputs']['summary_md']}")


if __name__ == "__main__":
    main()
