# Eval — bmv

- **Sample type:** BMV
- **Script:** `bmv_measure@1.0`
- **Reference:** `benchmarks/bmv/reference_runs.csv` (11 prior images)

> ## How to read this report
>
> This compares your run against the lab's reference set of past trusted runs.
> We check whether the **measurement quality** looks normal — *not* whether
> the absolute sizes are what you expected. Sizes are what the script
> measures; this report doesn't judge them. (Your VLP25 sample might turn out
> to be 20 nm. The eval doesn't care, as long as the measurement was clean.)
>
> ### What's a warning?
> Each per-image warning means a quality metric fell outside the typical range
> we see across past measurements. The "typical range" is **median ± 2 × IQR**
> (interquartile range — the spread of the middle 50% of past values).
> Roughly speaking, half the past runs sit within median ± 1 IQR; almost all
> within median ± 2 IQR. Outside that = unusual, worth a look.
>
> ### What to do when you see a warning
> Open the corresponding overlay PNG in `overlays/`. Most quality warnings
> show up visually — uneven staining, particle clustering, focus issues, weak
> contrast. The eye usually sees the cause within seconds. The numbers just
> tell you which image to look at.
>
> ### Metric reference

> - **`wall_fit_success_rate`** — Fraction of detected particles for which the script found a clean capsid wall. Low = the protein ring isn't clear enough in this image (stain too thin, contrast too low, or focus issues).
> - **`reliable_rate`** — Fraction of detected particles whose capsid measurement passed all quality filters (good wall fit + circular wall). Low = many particles look problematic to the algorithm.
> - **`drop_rate`** — 1 − reliable_rate. Fraction excluded from the final mean.
> - **`median_wall_cv`** — Per-particle, the script measures wall radius in 8 angular sectors and computes std/mean (= CV). 0 = perfect circle. The median CV across all particles in this image is reported. Higher = particles look non-circular (deformed walls, particle clustering, bad fits).
> - **`iqr_wall_cv`** — Spread of wall_cv across particles in this image. High IQR = some particles fit cleanly, others don't — often uneven staining.
> - **`capsid_mean_nm`** — Mean capsid diameter in nm across all reliable particles in this image. Gated only for BMV (where every prep is the same biological entity and drift in absolute size is a real signal). Not gated for VLP — different VLP samples have different sizes by design.
>
> Absolute capsid size **is** gated for this sample type (BMV is a single well-characterised virus species).

## Headline
- 11 images evaluated
- **0 per-image quality warning(s)**
- This run measured: capsid 28.59 ± 0.94 nm (reported and gated)
- Hand vs script: Δ capsid **+0.61 nm** (hand n=336 from 1 CSV(s))

## Per-image quality check

All images within reference range. ✓

## Hand vs script

Read **336 particles** from 1 CSV(s): bmv_hand_measurements.csv

| metric | hand | script | Δ (hand − script) |
|---|---|---|---|
| Capsid mean | 29.20 nm | 28.59 nm | **+0.61 nm** |

_A persistent positive or negative delta across runs would suggest the script reads systematically larger or smaller than human ImageJ tracing — a calibration issue rather than a per-batch problem._
