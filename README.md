# labflow-ai

TEM image-analysis co-scientist for a biology lab. Automates particle measurement (gold NPs, virus capsids, capsomeres) from Gatan `.dm3`/`.dm4` images, with results trustworthy enough to replace or augment hand measurement in ImageJ.

## What's here

| | |
|---|---|
| `vlp_measure_v2.py` | VLP gold + capsid measurement. BMV-style radial-profile wall fit, quality-based reliability, mag-invariant smoothing. **Primary VLP script.** |
| `vlp_measure.py` | Original VLP script (population-clip reliability). Kept for comparison; v2 is the one to use. |
| `bmv_measure.py` | BMV / BOG capsid measurement (no gold anchor). Hough detection + WALL_DESCENT_FRAC wall fit. |
| `plot_vlp_scatter.py` | Combined scatter / 2D-hist / KDE / pooled hist across multiple VLP samples. |
| `plot_capsid_groups.py` | Per-image strip plot showing gold + capsid medians side-by-side. |
| `LAB_NOTEBOOK.md` | Decision log: every non-trivial measurement decision, with dates and reasons. Read this first when revisiting. |
| `CLAUDE.md` | Co-scientist working principles (decompose, anchor on easiest feature, generate overlays, surface anomalies). |
| `.claude/skills/lab-pipeline/` | Project-scoped Claude Code skill. Documents the tool surface, file conventions, and standard scientist workflows. |
| `benchmarks/` | *(planned)* Gold-standard runs and hand measurements per sample type. Curated manually. |
| `incoming/` | *(convention)* DM3 dump location for new sample batches. Each batch in its own dated subfolder. |
| `results/` | Per-run outputs: CSVs, overlays, plots, `SUMMARY.md`. |

## How to run a measurement

Programmatic (preferred when an agent is driving):

```python
from vlp_measure_v2 import run

result = run(
    image_path="incoming/2026-05-15_VLP17_batch4",
    sample_type="VLP",
    out_dir="results/vlp17_batch4",
    workers=6,
)
print(result["summary"])  # n_reliable, gold/capsid mean ± std, drop rate, …
```

CLI (preferred when a human is driving):

```bash
uv run python vlp_measure_v2.py "incoming/2026-05-15_VLP17_batch4" \
    --sample-type VLP --out results/vlp17_batch4 --workers 6
```

Both paths produce the same outputs: `vlp_measurements.csv`, `overlays/*.png`, histograms, scatter, and a `SUMMARY.md` with the headline numbers + per-image table.

## Data model (current)

- `incoming/<date>_<sample>_batch<N>/` — DM3/DM4 dump, optional `sample.txt` describing sample type.
- `results/<run_name>/` — measurement outputs from one `run()`. Includes `SUMMARY.md` for the headline.
- `LAB_NOTEBOOK.md` — append-only log of substantive decisions.

Planned (designed in `LAB_NOTEBOOK.md`, not yet implemented):

- `benchmarks/<sample_type>/gold_standard.csv` — curated approved runs (manual gate).
- `benchmarks/<sample_type>/hand_measurements.csv` — per-particle hand data, joined to runs by `batch_id`.
- `eval.py` — script-vs-gold-script + script-vs-hand comparison engine.
- `approve.py` — manual gate for promoting a run into the gold-standard set.
- `validate_script.py` — re-run on the gold standard after a script change to catch regressions.

## Working principles

See `CLAUDE.md` for the full list. The non-negotiables:

1. **Always look at the actual image first** before picking a detector or diagnosing a bad result. Most mistakes come from theorising over CSV stats instead of opening the overlay.
2. **Generate overlay PNGs for every run** — numbers alone cannot validate a pipeline.
3. **Auto-threshold per image.** CLAHE and staining variation mean fixed cutoffs misbehave.
4. **Declare physical scales in nm**, convert to pixels per-image. Pixel-only constants cause magnification-dependent bias.
5. **Surface anomalies** — 100% detection rates, step changes, outliers in range. Print and flag, don't silently accept.

## Sample types

| Sample | Anchors on | Script |
|---|---|---|
| VLPs (gold-NP-cored) — VLP17, VLP20, VLP_100 | Gold NP (near-black, circular) | `vlp_measure_v2.py` |
| Bare gold NPs | Gold NP | `vlp_measure_v2.py` (gold path only) |
| BMV / BOG (no gold) | Bright protein ring, dark stain pool | `bmv_measure.py` |

Sample classification from a raw image is **not yet automated** — the scientist supplies it (via `sample.txt` or as a `run()` argument).

## Setup

```bash
uv sync           # install dependencies
uv run python vlp_measure_v2.py --help
```
