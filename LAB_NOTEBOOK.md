# Lab Notebook — labflow-ai

---

## 2026-04-03 — VLP TEM Image Analysis: Gold NP + Capsid Diameter

### Goal
Automate measurement of two quantities from negative-stain TEM images of virus-like particles (VLPs): the diameter of the encapsulated gold nanoparticle core, and the outer diameter of the surrounding capsid shell. Images are Gatan Digital Micrograph format (.dm3 / .dm4) with embedded pixel calibration.

### Sample
- **VLP17**: 17 images, gold NPs ~16-17 nm, imaged at 0.3536 nm/px (most) and 0.1429 nm/px (3 higher-magnification images)
- **VLP20**: 20 images, gold NPs ~20-21 nm, imaged at 0.3536 nm/px

### Approach

**Step 1: Gold NP detection**

Initial attempt used Difference-of-Gaussians (DoG) blob detection — ran for days on 4096×4096 images, abandoned. Switched to intensity thresholding + connected components, which ran in seconds.

Key issue: a fixed threshold failed across images because CLAHE redistributes contrast differently per image. Solved with per-image automatic threshold selection (Otsu on the dark portion of the histogram). This gave stable results across all images without manual tuning.

Detections filtered by circularity (solidity > 0.75, eccentricity < 0.75) and size (10–30 nm) to reject noise and aggregates.

**Step 2: Capsid detection**

Approach: radial density profiles (RDPs) sampled outward from each confirmed gold NP center.

Initial implementation searched for the steepest gradient (peak of the intensity rise). This landed consistently too far out — it was finding the outer edge of the capsid wall where it returns to background, not the capsid wall itself.

Diagnosis: saved RDP plots as images with vertical lines marking the detected gold and capsid radii. Particles were numbered in the overlay PNG to allow direct matching between overlay and profile. By identifying specific particles that looked correct in the overlay and examining their profiles, a clear pattern emerged: in correct detections, the capsid line sat at the **local minimum** of the smoothed profile at ~20 nm from center (for 17 nm gold NPs). The capsid wall is a dark protein ring — a minimum in the RDP, not a gradient peak.

Updated algorithm to search for the first significant local minimum using `scipy.signal.find_peaks` on the inverted smoothed profile, with:
- A 5 nm skip at the search start to avoid catching the intensity recovery immediately adjacent to the gold NP
- Prominence filtering (0.02) to reject shallow noise dips
- Sub-pixel parabolic interpolation around the minimum to avoid quantization artifacts

**Step 3: Quantization fix**

Scatter plot of capsid vs gold NP diameter showed horizontal banding — capsid measurements landing at discrete multiples of ~0.71 nm (2 × 0.354 nm/px pixel size). Sub-pixel parabolic interpolation around the profile minimum removed this.

### Results

| Sample | Gold NP mean ± std | Capsid mean ± std | n |
|--------|-------------------|-------------------|---|
| VLP17 | 16.8 ± 1.5 nm | 38.5 ± 2.7 nm | 4059 |
| VLP20 | 20.8 ± 1.3 nm | 43.4 ± 2.2 nm | 271 |

Detection rate: 100% of detected gold NPs received a capsid measurement, consistent with spot-check visual inspection showing all detected particles were encapsidated. Empty capsids (no gold NP) are not detected by this pipeline — detection is anchored on the gold NP.

A step change in capsid diameter between VLP17_0001–0006 (mean ~37 nm) and VLP17_0007–0013 (mean ~39 nm) was detected, likely reflecting different areas of the TEM grid imaged in the same session.

VLP17_0017 showed distinctly larger gold NPs (20.9 nm) — consistent with a different NP batch.

### Outputs per run
- `results/vlp_measurements.csv` — per-particle gold and capsid diameters
- `results/overlays/<name>_overlay.png` — detection overlay with numbered particle labels
- `results/vlp_histograms.png` — 1 nm bin histograms for gold and capsid
- `results/vlp_scatter.png` — scatter plot of capsid vs gold NP diameter
- `results/overlays/<name>_profiles.png` — radial profiles for selected particles (debug mode)

### Script
`vlp_measure.py` — no size parameters required. Gold NP size and capsid size discovered automatically.

```
uv run python vlp_measure.py "images/VLPs for machine learning project" --pattern "VLP17_*"
uv run python vlp_measure.py "images/VLPs for machine learning project/VLP17_0001.dm4" --debug --debug-indices 11 61 83
```

### Open questions
- Capsid detection for empty capsids (no gold NP anchor)
- Formal eval against ImageJ hand measurements
- Capsid detection for plain viruses (BMV) — no gold NP anchor
- ~~Whether the 103 nm capsid outlier is a real particle or a false detection~~ → resolved by reliability flagging

---

## 2026-04-08 — Cross-sample validation + Reliability Flagging

### Samples tested
- **VLP_100** (VLP_100_06_14_0001–0004.dm3): 1024×1024 px at 0.5865 nm/px — different pixel size and smaller gold NPs (~12 nm)

### Findings

Script ran on VLP_100 without parameter changes. Auto-threshold handled the smaller gold NPs and coarser pixel size correctly.

| Sample | Gold NP mean ± std | Capsid mean ± std (reliable) | n |
|--------|-------------------|------------------------------|---|
| VLP17 | 16.8 ± 1.5 nm | 38.5 ± 2.7 nm | 4059 |
| VLP20 | 20.8 ± 1.3 nm | 43.4 ± 2.2 nm | 271 |
| VLP_100 | 12.2 ± 0.8 nm | 34.5 ± 1.3 nm | 252 |

### Problem identified: false capsid detections

Spot-checking particles #81 and #90 in VLP_100_06_14_0002 showed oversized cyan circles in the overlay. RDP inspection revealed two distinct failure modes:

1. **Shallow dip** (#81): real capsid dip at ~15 nm was not prominent enough to be detected; algorithm found a weaker but more prominent feature further out. More likely to occur at coarser pixel sizes (0.5865 nm/px) where profiles are noisier.

2. **Contamination** (#90): a large dark aggregate adjacent to the particle contaminated the radial profile, creating a false deep minimum at ~31 nm radius.

### Fix: population-level reliability flagging

After measuring all particles per image, computed per-image median and std of capsid diameters. Particles with capsid diameter > 2 std from the median are flagged `is_reliable = False` in the CSV. Measurements are not modified — flagged particles are retained for transparency but excluded from summary stats and shown as red X marks in the scatter plot.

Effect on VLP_100_0002: 4 particles flagged, capsid std dropped from 3.9 nm → 1.3 nm, max capsid from 62 nm → 38 nm.

### Updated outputs
- `is_reliable` column added to CSV
- Scatter plot shows flagged particles as red X marks separately from reliable detections
- Summary prints reliable-only capsid stats alongside flag count

### Open questions
- Capsid detection for empty capsids (no gold NP anchor)
- Formal eval against ImageJ hand measurements
- Capsid detection for plain viruses (BMV) — no gold NP anchor
- Shallow dip failure mode: more smoothing at coarser pixel sizes, or lower prominence threshold with tighter search window?

---

## 2026-04-25 — BMV / BOG Capsid Measurement (no gold anchor)

### Goal
Measure outer capsid diameter for plain BMV and BOG (BMV labeled with a fluorescent protein) in negative-stain TEM, where there is no gold NP to anchor on. Target: 28 nm published BMV diameter.

### Sample
- **BMV_a–k**: 11 images at 0.0976, 0.234 nm/px (mixed magnifications)
- **BOG298 series**: 0.953 nm/px (much lower mag, denser packing)
- All uranyl-acetate negative-stained

### Key signature observation
BMV in this stain is **donut-shaped**: dark center (stain fills the interior) + bright ring (protein wall, stain excluded) + dark stain pool around the outside + lighter background. This is fundamentally different from VLP gold-anchored particles, where the gold NP is a near-black blob at the center.

### Approach evolution

**Initial attempt — LoG bright-blob detection**: failed catastrophically. LoG at the expected scale matches *bright* circular regions, but BMV has a *dark* center. Detector found bright stain texture and missed real BMVs entirely. Switched immediately.

**Hough circle transform**: works regardless of bright/dark interior — votes from intensity-gradient edges toward circle centers. With NMS at 1.4 × expected radius (= just-touching capsids), recovers densely-packed real particles cleanly.

**Center refinement**: Hough vote peaks can sit a few pixels off the true center, biasing the wall fit small. Snap each candidate to the local intensity minimum (dark stain core) within a small window — but **only keep the snap if it lowers per-sector wall_cv** (the snap can be misled into adjacent stain pools). This conditional-snap pattern was important for not destroying recall.

**Wall fit definition**: tested several:
- Steepest descent of radial profile (used in VLP code) — too small (~26 nm)
- Half-max between bright-ring peak and post-peak stain minimum — still too small visually
- Post-peak stain minimum directly — too big (~33 nm)
- **0.75 of the way from peak to post-peak min** — matches visually-traced outer edge (~28.6 nm overall mean)

The free parameter `WALL_DESCENT_FRAC = 0.75` controls this; 0.5 is half-max, 1.0 is the stain-minimum.

**Quality filters** — calibrated against user-labeled good (#349, #124) vs bad (#48, #153, #574, #313, #347) examples on BMV_a:

| filter | threshold | what it kills |
|---|---|---|
| `MIN_CONTRAST = 0.32` | ring_mean − exterior_mean | low-contrast aggregates / debris |
| `MAX_WALL_CV = 0.18` | per-sector wall radius coefficient of variation | non-circular stain texture (the strongest discriminator) |
| `MAX_EXTERIOR_MEAN = 0.32` | absolute exterior intensity | particles sitting in bright stain regions |
| `DIAM_TOL_FRAC = 0.30` | diameter range gate | fits gone wildly off-scale |

`MIN_UNIFORMITY` was tested but the good/bad classes overlapped — disabled.

**Overlap exclusion**: touching particles share a stain ring, so both walls converge on it. After filtering, pairs whose centres are closer than `(r1 + r2) × 0.95` are demoted to unreliable. Implemented with cKDTree pairwise queries.

**Magnification-invariant smoothing**: critical. Initial `PROFILE_SMOOTH_SIGMA = 3.0 px` produced a +1 nm bias at high mag (BMV_h–k at 0.0976 nm/px) vs low mag (BMV_a at 0.234 nm/px). Reason: 3 px = 0.29 nm at high mag vs 0.70 nm at low mag, so the bright-ring peak gets pulled inward more at low mag, shifting the wall measurement smaller. Fixed by declaring smoothing in nm (`PROFILE_SMOOTH_NM = 0.7`) and converting to px per-image. The 0.7 nm value is inherent to negative-stain virus analysis — it sits below the ~1.5 nm physical width of the protein-stain transition and above per-pixel noise scales.

### Visualisation discipline (re-affirmed)
- Cyan circles = reliable detections only by default; tomato + magenta hidden behind `--show-flagged`
- Numbered labels only on reliable particles, so debug profile plots match what's visible
- Per-image overlays mandatory; saved alongside CSV every run
- Per-image summary plot (`bmv_summary.png`) shows means ± std side-by-side with combined histogram so per-image drift can be spotted

### Result
**n = 1011 reliable across 11 BMV images. Mean 28.6 ± 1.0 nm. Range 24.8–31.8 nm.** After mag-invariant smoothing fix, the BMV_h–k high-mag bias should resolve to ~28 nm consistent with a–g.

### Parallelisation
Per-image work is independent — wrapped in `process_image(path, ...)` and dispatched via `concurrent.futures.ProcessPoolExecutor`. Default `--workers 4`. Smaller wall-clock by ~5× at 6 workers on the BMV set. Algorithm changes don't complicate parallelisation since the worker function takes only `(path, output dir, expected_nm, debug flags)` and returns `(per-image df, log string)`.

### Open questions / next steps
- Verify BMV_h–k bias is gone after mag-invariant smoothing (re-run pending)
- Run BOG298 set — stricter `MAX_EXTERIOR_MEAN` may need loosening for the brighter background regions there
- Capsomere measurement inside verified BMVs (12 pentamers + 20 hexamers, T=3) — bookmarked for after capsid measurement is fully trusted
- Hand-validate against ImageJ measurements on a subset

---

## 2026-05-03 — VLP measurement v2 + per-image quality diagnosis

### Goal
Port the BMV-style wall fit and quality gating into the VLP pipeline (`vlp_measure_v2.py`) and validate against the hand measurements collected in ImageJ.

### What changed v1 → v2
1. **Wall fit**: replaced `argmin(gradient)` (which sits at the inflection point inside the protein density) with a `WALL_DESCENT_FRAC = 0.75` interpolation between the bright protein-ring peak and the darker exterior. Matches what a human draws as the outer protein edge.
2. **Reliability**: replaced the v1 population clip (`within ±2σ of median capsid_diameter`) with per-particle quality gates: `ring_contrast ≥ 0.05`, `peak_prominence ≥ 0.15`, `wall_radius_cv ≤ 0.20`. These remove garbage fits whose numbers happen to land near the mean (which v1 silently kept).
3. **Magnification-invariant smoothing**: profile smoothing declared in nm (`PROFILE_SMOOTH_NM = 0.7`) and converted to px per-image — same fix as BMV.
4. **Overlap exclusion** via cKDTree pairwise distance, same as BMV.
5. **Plumbing**: `ProcessPoolExecutor` parallelism, tk progress GUI, `--show-flagged` overlay flag.

### Results (reliable particles only)
| sample | n | gold (nm) | capsid (nm) | hand capsid (nm) |
|---|---|---|---|---|
| VLP17  | 2964 | 16.7 ± 1.4 | 34.4 ± 1.8 | 34.4 |
| VLP20  | 1106 | 20.8 ± 1.3 | 39.0 ± 1.7 | – |
| VLP_100 | 288 | 12.1 ± 0.8 | 29.3 ± 1.0 | 28.4 |

VLP17 capsid mean lands within 0.05 nm of hand. VLP_100 reads ~1 nm bigger than hand — sample-specific tuning may be warranted but not pursued yet.

### Per-image quality breakdown (VLP17)
Per-image strip plot (`results/vlp17_v2/capsid_groups.png`) shows the gold NP median is rock-stable across all 17 images (~16.5–17.0 nm), but the capsid median drifts by ~1.5 nm across image groups. Drop-rate breakdown:

| group | images | gold | wall fit % | reliable % | drop % |
|---|---|---|---|---|---|
| 0001–0006 | 6 | 1257 | 67% | 54% | **46%** |
| 0007–0013 | 7 | 2509 | 95% | 83% | **17%** |
| 0014–0016 | 3 | 195 | 89% | 74% | 26% |
| 0017 | 1 | 99 | 76% | 66% | 34% (sample outlier — gold ~21 nm, capsid ~39 nm) |

### Diagnosis: stain uniformity, not staining intensity
Side-by-side overlays of VLP17_0003 (low yield, capsid 33.3 nm) vs VLP17_0010 (high yield, capsid 34.7 nm) show the same particle morphology — same dark gold cores, same bright protein rings — but the **background** differs:
- 0001–0006 fields show dark stain pooling into webs/cracks across the field. Each particle's radial profile has a noisy/biased baseline because the "outside-the-wall" intensity is itself splotchy, so the wall fit either fails (NaN) or comes out non-circular (high `wall_radius_cv`).
- 0007–0013 fields have a clean uniform background. The wall ring stands out cleanly against a flat baseline → fit nails it nearly every time.

This **also explains the ~1.5 nm capsid shift between groups** — but the direction is the *opposite* of what was first claimed (corrected below in the "Mechanism corrected" subsection). The brighter 0001–0006 baseline biases capsid measurements *small*; the darker 0007–0013 baseline lets the descent reach the true outer protein-stain boundary.

So both the drop rate AND the capsid size shift are downstream of the same artifact: stain uniformity at grid prep. Not an algorithm problem; not a magnification problem; not a biology problem.

### Mechanism corrected — direction of bias
Initial write-up had the staining direction flipped. Correct interpretation:

The wall fit reports the radius where intensity = `peak − WALL_DESCENT_FRAC × (peak − baseline)` = `0.25 × peak + 0.75 × baseline` (with `WALL_DESCENT_FRAC = 0.75`). The wall radius therefore depends on the local **baseline** (dark exterior) directly:
- **Darker baseline** → descent target sits at lower absolute intensity → descent walks **further outward** before triggering → **larger** measured radius.
- **Lighter baseline** → descent target sits closer to peak → descent stops **earlier** → **smaller** measured radius.

`bg_mean` per group: 0001–0006 = 0.533 (lighter); 0007–0013 = 0.509 (darker). Caveat: the σ = 100 nm low-pass doesn't fully mask particles, and 0007–0013 has more particles per field (~358 vs ~209), so part of the 0.024 `bg_mean` difference is particle density rather than stain. The wall fit itself reads a **per-particle local baseline** outside each particle's own wall, so it sees the local stain darkness directly — and that's the variable that's actually shifting the measurement.

So the corrected story:
- **0007–0013**: heavier, uniform uranyl acetate → dark, uniform exterior → descent reaches the true protein-stain boundary → **~34.6 nm (correct)**.
- **0001–0006**: lighter / patchier stain → brighter exterior on average → descent stops short → **~33.2 nm (under-measured)**.

This is also consistent with hand-measurement intuition — heavier negative stain gives crisper edges, easier for the eye, and the algorithm gets the right answer there too. Where stain is sparse, both eye and algorithm have less signal, but the algorithm's WALL_DESCENT_FRAC anchoring biases the result *systematically smaller*, not just noisier.

### Validation plan: hand-measure VLP17_0003 vs VLP17_0010
To confirm the bias is algorithmic (not real biology), measure capsids by hand in ImageJ on these two images:
- **VLP17_0003** (0001–0006 group): algorithm reports capsid median **33.3 nm** with high drop rate. Prediction: hand measurement should read closer to **34–35 nm**, matching the 0007–0013 algorithmic value. If hand confirms ~33 nm, the bias is real biology and the algorithm is right.
- **VLP17_0010** (0007–0013 group): algorithm reports **34.7 nm**, low drop rate, clean overlays. Prediction: hand should also read ~34.7 nm, confirming the algorithm's number is the unbiased one.

Outcome map:
| hand on 0003 | hand on 0010 | conclusion |
|---|---|---|
| ≈ 34.5 | ≈ 34.7 | Algorithm under-measures stain-poor fields. Need a stain-aware correction or per-image baseline normalisation. |
| ≈ 33.3 | ≈ 34.7 | Real ~1.5 nm size difference between grids — biology, not artifact. Unlikely given the bg_mean correlation but possible. |
| ≈ 33.3 | ≈ 33.3 | Algorithm has a systematic offset in 0007–0013 too — would need re-tuning. |

Save hand measurements to `results/vlp17/VLP17_0003_hand.csv` and `results/vlp17/VLP17_0010_hand.csv` (one row per particle, ImageJ Length export in nm). Plot hand vs algorithm per particle for both images.

### Hand-measurement results (validated)
Hand measurements collected in ImageJ on n=61 particles per image:

| image | hand mean | algo mean | algo − hand |
|---|---|---|---|
| VLP17_0003 (under-stained group) | 33.87 nm | 33.38 nm | **−0.49 nm** |
| VLP17_0010 (well-stained group)  | 34.64 nm | 34.78 nm | +0.14 nm |
| **Hand-measured gap between groups** | **0.77 nm** | algo gap: 1.40 nm | algo amplifies by ~0.6 nm |

Hand measurements files: `results/vlp17/vlp17_003_hand.csv`, `results/vlp17/vlp17_0010_hand.csv`.

### Conclusion
The inter-image capsid difference is downstream of **staining/grid-prep variation between imaging sessions**, not real biology and not an algorithm bug:
- ~0.8 nm of the gap is real (visible by eye and confirmed in the hand data) — likely staining depth or Fresnel-fringe defocus differences between grids.
- ~0.6 nm extra is algorithmic amplification: the wall fit anchors its descent target to the local baseline (`target = 0.25 × peak + 0.75 × baseline`), so when the local baseline differs between fields, the measured wall radius shifts in the same direction.
- The algorithm reads correctly (within 0.14 nm of hand) on cleanly-stained fields. It under-reads by ~0.5 nm on under-stained / locally-variable fields.

### Important caveat on the visual diagnosis
The PNG overlays we used to "see" stain uniformity were processed through `normalise()` = (per-image percentile stretch to [0,1]) + (CLAHE, `clip_limit=0.02`). Both stages amplify contrast, and CLAHE specifically boosts *local* contrast within tiles. So the dramatic "dark webbing in 0003 vs uniform grey in 0010" is partly a display artifact: tiny absolute differences in the raw .dm4 get exaggerated to look like big visual contrasts. Viewing the raw .dm4 in ImageJ shows much more modest underlying contrast in both images.

Implication: per-image bg metrics computed on the post-`normalise` image (the `bg_mean = 0.533 vs 0.509` numbers above) are **confounded by the per-image stretch** and shouldn't be read as "0001–0006 is brighter overall in real intensity terms." They describe the *shape* of each image's normalized distribution, not the absolute stain darkness. The `bg_std` and `bg_range` signals are more reliable than `bg_mean` because they describe spatial variability (which survives stretching) rather than absolute level (which doesn't).

The right next step here would be to compute the same diagnostics on raw 16-bit counts before `normalise()`, and possibly to run the wall fit on a globally-stretched (no-CLAHE) version to see if the inter-image bias shrinks. Tracked as a follow-up below.

### Quantitative support: per-image background metrics
To check the visual diagnosis against numbers, computed a low-pass-filtered background (Gaussian σ = 100 nm — much larger than any particle, so all particle structure is blurred away and only the slow background variation survives) and measured per-image:
- `bg_mean` — overall darkness of the field after particles are blurred out
- `bg_std`, `bg_range` (P95–P5) — how non-uniform the remaining background is

| group | bg mean | bg std | bg P95–P5 range | drop % | capsid (nm) |
|---|---|---|---|---|---|
| 0001–0006 | **0.533** | 0.020 | **0.066** | **46%** | 33.2 |
| 0007–0013 | 0.509 | 0.015 | 0.050 | 17% | 34.6 |
| 0014–0016 | 0.505 | 0.004 | 0.014 | 26% | 33.6 |
| 0017 | 0.515 | 0.015 | 0.048 | 34% | 39.1 |

Pearson correlations across the 17 images:
- `bg_mean` (overall darkness) vs drop rate: **r = 0.71**
- `bg_std` (unevenness) vs drop rate: r = 0.38
- `bg_range` (P95–P5) vs drop rate: r = 0.38

Both signals are present, but the dominant one is **overall darkness**: heavier-stained fields drop more particles. Unevenness is the secondary contributor. Together these support a single grid-prep story — "more stain deposited and it pooled non-uniformly" — rather than a pure brightness or pure unevenness story alone. The capsid-size shift between groups (0001–0006 reads ~33.2 nm vs 0007–0013 reads ~34.6 nm) tracks `bg_mean` in the same direction: darker exterior pulls the radial descent inward, biasing wall radius small.

(The 0014–0016 group reads tiny `bg_std` because it was acquired at higher magnification — 0.143 nm/px vs 0.354 — so the field of view is physically smaller and there's less long-wavelength background variation to capture at the 100-nm scale. The σ=100 nm low-pass is mag-invariant, but the field size relative to that filter scale is not.)

Plot saved to `results/vlp17_v2/stain_uniformity.png` (4-panel: drop rate vs `bg_std`, drop rate vs `bg_range`, drop rate vs `bg_mean`, capsid size vs `bg_range`).

### Implications for the data-analysis agent
- A per-image quality summary (`n_detected`, `n_reliable`, fraction passing each filter) is sufficient to flag uneven-stain images automatically — when wall-fit success drops below ~80% the per-image capsid mean should not be trusted absolutely (it's biased low).
- Staining unevenness shows up in the algorithm before it's visible to the eye on a single particle. The aggregate drop rate is a more sensitive uniformity indicator than visual inspection of one particle.

### Outputs
- `results/vlp17_v2/`, `results/vlp20_v2/`, `results/vlp100_v2/` — per-sample CSVs, histograms, scatter, overlays
- `results/combined_v2/` — combined scatter, 2D histogram, KDE, pooled capsid hist across all three samples
- `results/comparisons/` — hand vs v1 vs v2 distribution comparisons
- `results/archive_v1/` — old v1 loose files moved out of `results/` root

### Open questions / next steps
- ✅ Hand-measure VLP17_0003 and VLP17_0010 — done; conclusion: ~0.8 nm real grid-prep difference + ~0.6 nm algorithmic amplification on under-stained fields
- Recompute per-image stain diagnostics on **raw 16-bit counts** (pre-`normalise`) to remove the per-image stretch confound; see if `bg_mean` ordering reverses
- Try running the wall fit on a globally-stretched (no-CLAHE) image and check whether the inter-group capsid gap shrinks
- VLP_100 capsid +1 nm vs hand: tune `WALL_DESCENT_FRAC` per-sample, or accept as the cost of one universal value?
- Add per-image quality summary to CSV output as a first-class column for the agent to gate on (include `bg_mean`, `bg_std`, `bg_range` from the σ = 100 nm low-pass)
- Define the "minimum quality bar" thresholds (e.g. wall-fit success ≥ 80%, n_reliable ≥ 30) per sample type

---

## 2026-05-08 — BMV speedup options (future work)

After getting `validate_script.py` running on the lean reference set, BMV is dominating runtime: ~75–115 s/image locally on M-series, ~2× that on x86 CI runners. VLP is 9–22 s/image (or 1–3 s for the smaller VLP_100 set). Hough circle transform on 4096×4096 images is the bottleneck.

Ranked by bang-for-buck if/when speedup becomes a priority:

### Quick wins (a few hours each)
1. **Downsample for Hough detection (3–5× BMV speedup).** Hough scales as O(W·H·N_radii). Find candidate centers on a 1024×1024 (4× downsampled) image, then refine the wall fit at full 4k. We need 4k for measuring the wall edge precisely, not for finding particle centers. Decoupling the two is the cleanest single optimization. Validate with `validate_script.py` to confirm same particles are found.
2. **Switch `skimage.transform.hough_circle` → `cv2.HoughCircles`** (5–10× on the same workload; OpenCV's C++ impl). Roughly drop-in; need to verify candidate set matches.
3. **Vectorize per-sector wall fitting** (2–3× on that step, applies to both VLP and BMV). Currently 8 angular sectors × ~30 radii × 360 angle bins via Python `for` loops; replace with NumPy.
4. **Skip overlay generation when `validate_script` is calling `process_image`.** Already an option in the design — just need to thread a `save_overlays=False` flag through. Saves ~1 s/image and avoids matplotlib in the hot path entirely.

### Medium projects (1–2 days)
5. **Replace CLAHE with simpler percentile stretch for BMV.** CLAHE matters for VLPs (gold contrast). For BMV, the protein–stain contrast is naturally high; a percentile stretch may be sufficient. Saves 5–10 s/image. Risk: quality filters are tuned to CLAHE-equalised images; would need re-validation.
6. **Reuse worker processes across multiple images.** Currently each parallel worker spawns, imports `bmv_measure` (matplotlib + scipy + skimage), processes one image, exits — repeating the ~1–2 s import cost per image. Persistent workers via `multiprocessing.Pool.imap` with `chunksize > 1` cut that overhead.

### Bigger projects (only if real volume demands it)
7. **GPU acceleration** for Hough + CLAHE via `cv2.cuda.HoughCircles` or a small PyTorch port. 10–50× on CUDA hardware. Adds a GPU dependency; right answer only if the lab moves to cloud / Modal-style execution.
8. **Replace Hough with a learned detector** (small U-Net or YOLO trained on existing BMV reference images). Faster *and* potentially more robust to weird stain conditions. Bigger lift; only worth it if BMV throughput becomes a real bottleneck.

### Recommended order if/when prioritised
Do **#1** first — biggest single win, conceptually clean (separation of "where" vs "exactly how big"), and `validate_script.py` will tell us immediately if it shifted measurements. With #1 alone, BMV likely drops from 75 s/image to 15–20 s/image; the lean reference set goes from ~4 min → ~1 min locally and ~6 min → ~2 min on CI.

**#3** (vectorised wall fitting) is the next best because it benefits both VLP and BMV and is purely a refactor.

Don't pursue **#7** or **#8** until BMV throughput is a real bottleneck. ROI on the downsampled Hough is much higher than on infra migration.

---

## 2026-05-08 — Validation made manual (CI cost)

Auto-running `validate.yml` on every PR cost ~10 min of CI per run. Most PRs don't change measurement behaviour, so this was wasted compute.

Switched the `validate.yml` workflow to **`workflow_dispatch` only** — manual trigger from the Actions tab. Added a separate lightweight `validate-reminder.yml` that posts a sticky comment on PRs touching `analysis/`, `benchmarks/`, or deps, reminding the author to either:

- Run validate locally (~3–4 min on M-series), or
- Trigger the manual workflow on the Actions tab.

Net effect: zero auto-CI minutes spent on regression testing; reviewers still see a visible reminder when it's worth running.

If BMV speedup work (above) lands and the regression run drops to ~2 min on CI, revisit and consider re-enabling auto-trigger on `pull_request` for `analysis/**` paths.
