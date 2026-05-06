# labflow-ai

TEM image-analysis co-scientist for a biology lab. Automates particle measurement (gold NPs, virus capsids, capsomeres) from Gatan `.dm3`/`.dm4` images, with results trustworthy enough to replace or augment hand measurement in ImageJ.

The whole stack works equally well via plain CLI (a scientist running scripts directly) and via Claude / Codex (an agent driving the same scripts conversationally). Claude is a convenience, never a dependency.

## What's here

| | |
|---|---|
| `analysis/vlp_measure_v2.py` | VLP gold + capsid measurement. **Self-contained — drop just this file alongside images and run it.** Image loading, normalisation, gold detection, capsid wall fit, overlays, summary. |
| `analysis/bmv_measure.py` | BMV / BOG capsid measurement (no gold anchor). Hough detection + WALL_DESCENT_FRAC wall fit. |
| `analysis/plot_vlp_scatter.py` | Combined scatter / 2D-hist / KDE / pooled hist across multiple VLP samples. |
| `analysis/plot_capsid_groups.py` | Per-image strip plot showing gold + capsid medians side-by-side. |
| `analysis/eval.py` | Compares a measurement run against the reference benchmarks (per-image quality + hand vs script). |
| `analysis/seed_benchmarks.py` | Initial seeder for `benchmarks/<sample_type>/`. |
| `analysis/add_to_reference.py` | Manual gate for appending new approved runs / hand data into the benchmarks. |
| `benchmarks/<sample_type>/` | `reference_runs.csv` (script outputs we trust) + `reference_hand.csv` (per-particle hand measurements). Eval compares new runs against these. |
| `LAB_NOTEBOOK.md` | Decision log: every non-trivial measurement decision, with dates and reasons. Read this first when revisiting. |
| `CLAUDE.md` / `AGENTS.md` | Co-scientist working principles + auto-bootstrap rule for Claude / Codex. |
| `.claude/skills/lab-pipeline/` | Project-scoped Claude Code skill. Tool surface, file conventions, scientist workflows. |
| `bootstrap.sh` | One-shot installer for new scientist repos. |
| `incoming/` | *(convention)* DM3 dump location for new sample batches. |
| `results/` | Per-run outputs: CSVs, overlays, plots, `SUMMARY.md`, `eval_report.md`. |

## Three install paths (pick whichever matches how you work)

### 1. Agent-driven install — the easiest path, no shell knowledge required

If you use Claude Code or Codex (CLI or desktop), you don't need to know `curl` or `bash`. Just give the agent the one file it reads at session start:

1. Make a fresh empty folder where you want to do TEM analysis.
2. Save the agent context file into it. Open https://github.com/lilydelalande/labflow-ai/blob/main/CLAUDE.md in your browser, click *Raw*, and save into your folder as:
   - `CLAUDE.md` — if you use Claude Code
   - `AGENTS.md` — if you use Codex (it's the same file content, just under a different name)

   (Or, from the terminal, `curl -sSL https://raw.githubusercontent.com/lilydelalande/labflow-ai/main/CLAUDE.md -o <CLAUDE.md or AGENTS.md>`.)
3. Open Claude Code or Codex in that folder and say *hi*.
4. The agent reads the context file, notices the analysis tools aren't installed, asks for your permission, and runs the bootstrap itself.

After that, drop your DM3s into `incoming/<batch_name>/` and ask the agent to analyze them.

### 2. Manual install — for technical users who prefer doing it themselves

```bash
cd ~/my-tem-project
curl -sSL https://raw.githubusercontent.com/lilydelalande/labflow-ai/main/bootstrap.sh | sh
```

This clones the repo into a hidden `.labflow/` cache, symlinks `analysis/`, `benchmarks/`, and the lab-pipeline skill into the working directory, copies `CLAUDE.md` (and symlinks `AGENTS.md` to it), gitignores the cache + data folders, runs `uv sync` to install Python deps, and creates empty `incoming/` and `results/` folders.

To update later: re-run the same one-liner — it's safe to re-run.

### 3. Single-file usage — just the measurement script, no agent or benchmarks

`analysis/vlp_measure_v2.py` is self-contained. If you only want VLP measurement and don't need eval / benchmarks / agent layers, drop the file alongside a folder of images and run:

```bash
mkdir my-tem-project && cd my-tem-project
curl -sSL https://raw.githubusercontent.com/lilydelalande/labflow-ai/main/analysis/vlp_measure_v2.py -o vlp_measure_v2.py
uv init && uv add ncempy pandas matplotlib scipy scikit-image tqdm
mv ~/my-images.zip . && unzip my-images.zip   # or however your DM3s arrive
uv run python vlp_measure_v2.py my-images/ --workers 6
```

Outputs land in `results/vlp_v2/` (the script creates the folder if it doesn't exist; you don't need to make it yourself). The eval auto-import skips silently when the broader stack isn't around.

---

After any of paths 1 or 2, the scientist's daily flow is the same: drop DM3s into `incoming/<batch>/`, then either:
- **Talk to Claude / Codex:** "analyze the new images" — the lab-pipeline skill takes over
- **Run scripts directly from terminal:** `uv run python -m analysis.vlp_measure_v2 incoming/<batch> --workers 6`

Both produce the same outputs: CSV, overlays, plots, `SUMMARY.md`, and `eval_report.md`.

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

Both produce, inside the directory passed as `out_dir` (or `--out`):

```
results/<batch_name>/
├── SUMMARY.md            # headline + per-image table — read this first
├── vlp_measurements.csv  # one row per detected gold NP, with per-particle metrics
├── vlp_histograms.png    # gold + capsid distributions
├── vlp_scatter.png       # gold-vs-capsid scatter
└── overlays/
    ├── VLP17_0001_overlay.png
    └── …                 # one per input image, raw on left, detections on right
```

All paths in the dict returned by `run()` are absolute and point inside `out_dir`. By convention, mirror the incoming folder name: `incoming/foo_batch/` → `out_dir=results/foo_batch/`. `out_dir` is created if it doesn't exist.

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
