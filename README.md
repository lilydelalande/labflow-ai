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
| `CLAUDE.md` | Co-scientist working principles + auto-bootstrap rule. Codex users: rename to `AGENTS.md` once after install (Codex reads only that name). |
| `.claude/skills/lab-pipeline/` | Project-scoped Claude Code skill. Tool surface, file conventions, scientist workflows. |
| `bootstrap.sh` | One-shot installer for new scientist repos. |
| `incoming/` | *(convention)* DM3 dump location for new sample batches. |
| `results/` | Per-run outputs: CSVs, overlays, plots, `SUMMARY.md`, `eval_report.md`. |

## Three install paths (pick whichever matches how you work)

### 1. Agent-driven install — the easiest path, no shell knowledge required

If you use Claude Code or Codex (CLI or desktop), you don't need to know `curl` or `bash`. Just give the agent the one file it reads at session start:

1. Make a fresh empty folder where you want to do TEM analysis.
2. Download `CLAUDE.md` into it. From the terminal:
   ```bash
   curl -sSL https://raw.githubusercontent.com/lilydelalande/labflow-ai/main/CLAUDE.md -o CLAUDE.md
   ```
   Or in a browser: open https://github.com/lilydelalande/labflow-ai/blob/main/CLAUDE.md, click *Raw*, save as `CLAUDE.md` into your folder.
3. **If you use Codex**, rename it: `mv CLAUDE.md AGENTS.md`. Codex reads only `AGENTS.md`. Claude Code users skip this step.
4. Open Claude Code or Codex in that folder and say *hi*.
5. The agent reads the context file, notices the analysis tools aren't installed, asks for your permission, and runs the bootstrap itself.

After that, drop your DM3s into `incoming/<batch_name>/` and ask the agent to analyze them.

### 2. Manual install — for technical users who prefer doing it themselves

```bash
cd ~/my-tem-project
curl -sSL https://raw.githubusercontent.com/lilydelalande/labflow-ai/main/bootstrap.sh | sh
```

This clones the repo into a hidden `.labflow/` cache, symlinks `analysis/`, `benchmarks/`, and the lab-pipeline skill into the working directory, copies `CLAUDE.md`, gitignores the cache + data folders, runs `uv sync` to install Python deps, and creates empty `incoming/` and `results/` folders. (Codex users: `mv CLAUDE.md AGENTS.md` once after install.)

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

## CLI reference (all scripts, all flags)

Every script in `analysis/` is invocable directly. Canonical form is `uv run python -m analysis.<name>` (works from any cwd as long as you're in a labflow-ai repo). Direct-path form `uv run python /path/to/analysis/<name>.py` also works for the measurement scripts.

### Measurement scripts

#### `analysis.vlp_measure_v2` — VLP gold + capsid measurement

```bash
uv run python -m analysis.vlp_measure_v2 <image_path> [flags]
```

| flag | default | what it does |
|---|---|---|
| `image_path` (positional) | required | A `.dm3`/`.dm4` file or a folder of them |
| `--pattern` | `*` | Glob when `image_path` is a folder |
| `--out` | `results/vlp_v2` | Output directory (`SUMMARY.md`, CSV, overlays, etc. land here) |
| `--sample-type` | `VLP` | Sample type for routing eval — leave at default for VLPs |
| `--sample-name` | folder name | Override the auto-derived sample name |
| `--workers` | `4` | Parallel worker processes |
| `--gold-threshold` | auto (Otsu) | Override the per-image gold-detection intensity cutoff |
| `--min-gold-nm` | `7.0` | Lower size cutoff for gold detection |
| `--max-gold-nm` | `30.0` | Upper size cutoff for gold detection |
| `--show-flagged` | off | Also draw flagged-unreliable circles in the overlay (tomato) |
| `--gui` | off | Pop a Tk progress window (closing it doesn't stop the run) |

#### `analysis.bmv_measure` — BMV / BOG capsid measurement (no gold anchor)

```bash
uv run python -m analysis.bmv_measure <image_path> [flags]
```

| flag | default | what it does |
|---|---|---|
| `image_path` (positional) | required | A `.dm3`/`.dm4` file or a folder of them |
| `--pattern` | `*` | Glob when `image_path` is a folder |
| `--out` | `results/bmv` | Output directory |
| `--expected-nm` | `28.0` | Expected capsid diameter in nm — Hough searches around this |
| `--workers` | `4` | Parallel worker processes |
| `--show-flagged` | off | Draw flagged (tomato) and no-fit (magenta) circles in overlay |
| `--debug` | off | Save radial profile plots for a sample of particles |
| `--debug-indices N1 N2 …` | random sample | Specific particle indices to plot in debug mode |
| `--gui` | off | Pop a Tk progress window |

### Eval

#### `analysis.eval` — compare a run to the reference benchmarks

```bash
uv run python -m analysis.eval <run_dir> [flags]
```

| flag | default | what it does |
|---|---|---|
| `run_dir` (positional) | required | Directory with the measurement CSV (e.g. `results/<sample>`) |
| `--sample-type` | `VLP` | `VLP` or `BMV` — decides which benchmarks dir + which metrics get gated |
| `--benchmarks-dir` | `benchmarks/` | Override location of reference CSVs |
| `--no-report` | off | Skip writing `eval_report.md`; still returns the eval dict |
| `--json` | off | Print full result as JSON instead of the headline |

### Reference-set management (manual gate)

#### `analysis.seed_benchmarks` — one-shot initial seed

```bash
uv run python -m analysis.seed_benchmarks
```

No flags. Reads `results/vlp17_v2/`, `results/vlp20_v2/`, `results/vlp100_v2/`, and `results/bmv/`, writes `benchmarks/<sample_type>/{reference_runs,reference_hand}.csv`. Re-running overwrites — use `add_to_reference` for incremental growth instead.

#### `analysis.add_to_reference run` — append a run to `reference_runs.csv`

```bash
uv run python -m analysis.add_to_reference run <run_dir> --approver <name> [flags]
```

| flag | default | what it does |
|---|---|---|
| `run_dir` (positional) | required | Directory with `vlp_measurements.csv` or `bmv_measurements.csv` |
| `--sample-type` | `VLP` | `VLP` or `BMV` |
| `--sample-name` | folder name | Override the auto-derived sample name |
| `--approver` | required | Person approving this run (audit trail) |
| `--notes` | empty | Free-form note explaining why this run is reference-worthy |
| `--force` | off | Replace existing rows on `(sample_name, filename)` collision instead of refusing |
| `--benchmarks-dir` | `benchmarks/` | Override target directory |

#### `analysis.add_to_reference hand` — append hand measurements to `reference_hand.csv`

```bash
uv run python -m analysis.add_to_reference hand <hand_csv> --sample-name <name> --scientist <name> [flags]
```

| flag | default | what it does |
|---|---|---|
| `hand_csv` (positional) | required | ImageJ-exported CSV (paired gold+capsid for VLP, single column for BMV) |
| `--sample-name` | required | Identifier shared with the script run this hand data validates |
| `--sample-type` | `VLP` | `VLP` or `BMV` |
| `--unit` | `um` | Length unit in source CSV: `um` or `nm` |
| `--scientist` | required | Person who hand-measured (audit trail) |
| `--measure-date` | today | ISO date when measurements were taken |
| `--notes` | empty | Free-form note |
| `--run-dir` | none | Also copy the hand CSV into `<run-dir>/hand/` so per-batch eval picks it up |
| `--force` | off | Replace existing rows on `(sample_name, source_file)` collision |
| `--benchmarks-dir` | `benchmarks/` | Override target directory |

### Plotting helpers

#### `analysis.plot_vlp_scatter` — combined scatter / 2D-hist / KDE / pooled hist

```bash
uv run python -m analysis.plot_vlp_scatter <csv> [<csv> …] [flags]
```

| flag | default | what it does |
|---|---|---|
| `csvs` (positional, ≥1) | required | One or more `vlp_measurements.csv` files |
| `--out-dir` | first CSV's folder | Where the combined plots land |
| `--reliable-only` | on | Plot only reliable detections (default behavior) |
| `--combined-hist2d` | off | One pooled 2D histogram instead of per-sample subplots |

#### `analysis.plot_capsid_groups` — per-image strip plot, gold + capsid medians side-by-side

```bash
uv run python -m analysis.plot_capsid_groups <csv> [flags]
```

| flag | default | what it does |
|---|---|---|
| `csv` (positional) | required | A `vlp_measurements.csv` |
| `--out` | CSV's folder | Where `capsid_groups.png` lands |

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
