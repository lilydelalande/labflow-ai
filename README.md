# labflow-ai

TEM image-analysis co-scientist for a biology lab. Automates particle measurement (gold NPs, virus capsids, capsomeres) from Gatan `.dm3`/`.dm4` images, with results trustworthy enough to replace or augment hand measurement in ImageJ.

The whole stack works equally well via plain CLI (a scientist running scripts directly) and via Claude / Codex (an agent driving the same scripts conversationally). Claude is a convenience, never a dependency.

## What's here

| | |
|---|---|
| `analysis/vlp_measure_v2.py` | VLP gold + capsid measurement. BMV-style radial-profile wall fit, quality-based reliability, mag-invariant smoothing. **Primary VLP script.** |
| `analysis/vlp_measure.py` | Original VLP script (population-clip reliability). Kept for comparison; v2 is the one to use. |
| `analysis/bmv_measure.py` | BMV / BOG capsid measurement (no gold anchor). Hough detection + WALL_DESCENT_FRAC wall fit. |
| `analysis/plot_vlp_scatter.py` | Combined scatter / 2D-hist / KDE / pooled hist across multiple VLP samples. |
| `analysis/plot_capsid_groups.py` | Per-image strip plot showing gold + capsid medians side-by-side. |
| `LAB_NOTEBOOK.md` | Decision log: every non-trivial measurement decision, with dates and reasons. Read this first when revisiting. |
| `CLAUDE.md` | Co-scientist working principles + the auto-bootstrap rule any agent should follow. |
| `.claude/skills/lab-pipeline/` | Project-scoped Claude Code skill. Tool surface, file conventions, scientist workflows. |
| `bootstrap.sh` | One-shot installer for new scientist repos. |
| `benchmarks/` | *(planned)* Gold-standard runs and hand measurements per sample type. Curated manually. |
| `incoming/` | *(convention)* DM3 dump location for new sample batches. |
| `results/` | Per-run outputs: CSVs, overlays, plots, `SUMMARY.md`. |

## Installing the stack into a new lab repo

For a scientist setting up TEM analysis in a fresh directory:

```bash
cd ~/my-tem-project
curl -sSL https://raw.githubusercontent.com/lily-de/labflow-ai/main/bootstrap.sh | sh
```

This:
1. Clones labflow-ai into a hidden `.labflow/` cache (gitignored)
2. Symlinks `analysis/`, `benchmarks/`, and the lab-pipeline skill into the working directory
3. Copies `CLAUDE.md` (so the scientist can edit it locally)
4. Pins the install to the current upstream SHA in `.labflow/INSTALLED_SHA`

After that, the scientist can drop DM3s into `incoming/<batch>/` and either:
- **Talk to Claude / Codex:** "analyze the new images" — the lab-pipeline skill takes over
- **Run scripts directly from terminal:** `uv run python -m analysis.vlp_measure_v2 incoming/<batch> --workers 6`

Both paths produce the same outputs: CSV, overlays, plots, `SUMMARY.md`.

To update later: `cd .labflow && git pull && cd .. && .labflow/bootstrap.sh --relink` (or just re-run the curl one-liner — it's idempotent).

## Running a measurement (programmatic vs CLI)

**Programmatic — preferred when an agent is driving:**

```python
from analysis.vlp_measure_v2 import run

result = run(
    image_path="incoming/2026-05-15_VLP17_batch4",
    sample_type="VLP",
    out_dir="results/vlp17_batch4",
    workers=6,
)
print(result["summary"])  # n_reliable, gold/capsid mean ± std, drop rate, …
```

**CLI — preferred when a human is driving:**

```bash
uv run python -m analysis.vlp_measure_v2 "incoming/2026-05-15_VLP17_batch4" \
    --sample-type VLP --out results/vlp17_batch4 --workers 6
```

Both produce: `vlp_measurements.csv`, `overlays/*.png`, histograms, scatter, and a `SUMMARY.md` with the headline + per-image table.

## Data model

- `incoming/<batch_name>/` — DM3/DM4 dump, one folder per imaging session/grid.
- `results/<batch_name>/` — measurement outputs. Mirrors the incoming folder name 1:1.
- `LAB_NOTEBOOK.md` — append-only log of substantive decisions.

Planned (designed in `LAB_NOTEBOOK.md`, not yet implemented):

- `benchmarks/<sample_type>/gold_standard.csv` — curated approved runs (manual gate).
- `benchmarks/<sample_type>/hand_measurements.csv` — per-particle hand data, joined to runs by `batch_id`.
- `analysis/eval.py` — script-vs-gold-script + script-vs-hand comparison engine.
- `analysis/approve.py` — manual gate for promoting a run into the gold-standard set.
- `analysis/validate_script.py` — re-run on the gold standard after a script change to catch regressions.

## Working principles

See `CLAUDE.md` for the full list. Non-negotiables:

1. **Always look at the actual image first** before picking a detector or diagnosing a bad result.
2. **Generate overlay PNGs for every run** — numbers alone cannot validate a pipeline.
3. **Auto-threshold per image.** CLAHE and staining variation mean fixed cutoffs misbehave.
4. **Declare physical scales in nm**, convert to pixels per-image.
5. **Surface anomalies** — 100% detection rates, step changes, outliers. Print and flag.

## Sample types

| Sample | Anchors on | Script |
|---|---|---|
| VLPs (gold-NP-cored) — VLP17, VLP20, VLP_100 | Gold NP (near-black, circular) | `analysis/vlp_measure_v2.py` |
| Bare gold NPs | Gold NP | `analysis/vlp_measure_v2.py` (gold path only) |
| BMV / BOG (no gold) | Bright protein ring, dark stain pool | `analysis/bmv_measure.py` |

Sample classification from a raw image is **not yet automated** — the scientist supplies it.

## Setup (development on this repo itself)

```bash
uv sync                                                   # install dependencies
uv run python -m analysis.vlp_measure_v2 --help          # verify install
```
