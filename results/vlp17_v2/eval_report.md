# Eval — vlp17_v2

- **Sample type:** VLP
- **Script:** `vlp_measure_v2@2.0(wall=0.75;contrast=0.05;cv=0.2;prom=0.15;smooth=0.7;end=25.0)`
- **Reference:** `benchmarks/vlp/reference_runs.csv` (41 prior images)

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

> - **`wall_fit_success_rate`** — Fraction of detected gold particles for which the script found a clean capsid wall. Low = the protein ring isn't clear enough in this image (stain too thin, contrast too low, or focus issues).
> - **`reliable_rate`** — Fraction of detected gold particles whose capsid measurement passed all quality filters (good wall fit + circular wall). Low = many particles look problematic to the algorithm.
> - **`drop_rate`** — 1 − reliable_rate. Fraction excluded from the final mean.
> - **`median_wall_cv`** — Per-particle, the script measures wall radius in 8 angular sectors and computes std/mean (= CV). 0 = perfect circle. The median CV across all particles in this image is reported. Higher = particles look non-circular (deformed walls, particle clustering, bad fits).
> - **`iqr_wall_cv`** — Spread of wall_cv across particles in this image. High IQR = some particles fit cleanly, others don't — often uneven staining.
>
> Absolute sizes (`capsid_mean_nm`, `gold_mean_nm`, etc.) are reported but **not gated**.

## Headline
- 17 images evaluated
- **1 per-image quality warning(s)**
- This run measured: gold 16.70 ± 1.39 nm, capsid 34.42 ± 1.85 nm (reported, not judged)
- Hand vs script: no hand CSVs in this run's `hand/` folder

## Per-image quality check

| filename | warnings |
|---|---|
| `VLP17_0004.dm4` | `iqr_wall_cv` = 0.146 (ref median 0.045, IQR 0.042; Δ +0.101) |

_Reference distribution built from 41 prior images of this sample type._

## Hand vs script

No hand-measurement CSVs found in `<run_dir>/hand/`. Drop one or more ImageJ CSVs there and re-run eval to see the script-vs-hand comparison.
