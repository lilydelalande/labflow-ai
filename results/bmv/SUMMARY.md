# Measurement run — BMV

- **Sample:** `bmv`
- **Script:** `bmv_measure@1.0`
- **Run dir:** `results/bmv`

## Headline
- 11 images, 3986 particle detections
- **1456 reliable** (36.5%) — 2530 dropped (63.5%)
- **Capsid:** mean 28.59 ± 0.94 nm (median 28.60)

## Eval
- Per-image quality check: compared against 11 reference runs; **0 warning(s)**
- Hand vs script: Δ capsid **+0.61 nm** (hand n=336)
- Full report: `results/bmv/eval_report.md`

## Outputs
- CSV: `results/bmv/bmv_measurements.csv`
- Overlays: `results/bmv/overlays`
- Histograms: `results/bmv/bmv_histogram.png`

## Per-image
| filename | n_detections | n_reliable | reliable_rate | wall_fit_success_rate | capsid median (nm) |
|---|---|---|---|---|---|
| BMV_a.dm3 | 606 | 181 | 30% | 75% | 28.81 |
| BMV_b.dm3 | 520 | 247 | 48% | 83% | 28.31 |
| BMV_c.dm3 | 449 | 239 | 53% | 83% | 28.62 |
| BMV_d.dm3 | 535 | 180 | 34% | 82% | 29.01 |
| BMV_e.dm3 | 466 | 174 | 37% | 79% | 29.62 |
| BMV_f.dm3 | 506 | 170 | 34% | 90% | 27.86 |
| BMV_g.dm3 | 483 | 185 | 38% | 86% | 28.00 |
| BMV_h.dm3 | 110 | 24 | 22% | 72% | 29.57 |
| BMV_i.dm3 | 110 | 16 | 15% | 68% | 29.27 |
| BMV_j.dm3 | 110 | 16 | 15% | 81% | 29.92 |
| BMV_k.dm3 | 91 | 24 | 26% | 73% | 30.04 |

## Glossary

- **`reliable_rate`** — percentage of detected particles whose capsid measurement passed all quality filters. Low = many particles look problematic.
- **`wall_fit_success_rate`** — percentage of detected particles for which the script successfully fit a capsid wall. Low = the protein ring isn't clear enough (weak stain, poor contrast, focus issues).
- **CV (coefficient of variation)** — `std / mean`, a dimensionless way to express spread. Used in the eval as `median_wall_cv` (per particle, std/mean of wall radius across 8 angular sectors — 0 = perfect circle, higher = non-round wall) and `iqr_wall_cv` (spread of those CVs across particles in the image — high = some clean fits, some bad).
- See `results/bmv/eval_report.md` for full metric definitions and how warnings are flagged.
