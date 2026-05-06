# labflow-ai

TEM image-analysis co-scientist for a biology lab. Automates particle measurement (gold NPs, virus capsids, capsomeres) from Gatan `.dm3`/`.dm4` images, with results trustworthy enough to replace or augment hand measurement in ImageJ.

The same stack works two ways: **drive it with Claude or Codex** (conversational), or **run the scripts yourself** (terminal). Pick whichever fits how you work.

## Get started

### With Claude or Codex

1. Make a fresh empty folder.
2. Download `CLAUDE.md` into it:
   ```bash
   curl -sSL https://raw.githubusercontent.com/lilydelalande/labflow-ai/main/CLAUDE.md -o CLAUDE.md
   ```
   *Codex users:* `mv CLAUDE.md AGENTS.md` (Codex reads only `AGENTS.md`).
3. Open Claude Code or Codex in that folder and say *hi*.
4. The agent will detect the stack isn't installed yet, ask for permission, and run the bootstrap. Approve it.

### On your own (no agent)

Two options. Both end in the same place.

**Option A — clone the repo:**
```bash
git clone https://github.com/lilydelalande/labflow-ai
cd labflow-ai
uv sync
```

**Option B — run the bootstrap script** (same thing the agent runs, into a working directory of your choice):
```bash
mkdir my-tem-project && cd my-tem-project
curl -sSL https://raw.githubusercontent.com/lilydelalande/labflow-ai/main/bootstrap.sh | sh
```

Both give you the same setup: `analysis/`, `benchmarks/`, the lab-pipeline skill, dependencies installed, `incoming/` and `results/` ready for data. To update later: `git pull` (Option A) or re-run the curl one-liner (Option B). Both are idempotent.

## Daily flow

Drop DM3 / DM4 files into `incoming/<sample_name>/` (the folder name becomes the run identifier). Then either:

- **Talk to the agent:** "analyze the new images in incoming/foo" — the lab-pipeline skill takes over.
- **Run a script yourself:** `uv run python -m analysis.vlp_measure_v2 incoming/foo --workers 6` (or `analysis.bmv_measure` for BMV).

Outputs land in `results/<sample_name>/`. By convention input and output folder names match.

## What a run produces

```
results/<sample_name>/
├── SUMMARY.md            # headline + per-image table — read this first
├── eval_report.md        # quality + hand-vs-script (when benchmarks/<sample_type>/ exists)
├── vlp_measurements.csv  # one row per detected particle  (BMV: bmv_measurements.csv)
├── vlp_histograms.png    # gold + capsid distributions    (BMV: bmv_histogram.png)
├── vlp_scatter.png       # gold-vs-capsid scatter         (BMV: omitted)
└── overlays/             # per-image PNG; raw on left, detections on right
```

To enable the script-vs-hand comparison in `eval_report.md`, drop your ImageJ-exported CSV into `results/<sample_name>/hand/` (any filename) — eval reads it on the next run.

## Sample types

| `sample_type` | example samples | anchors on | script |
|---|---|---|---|
| `VLP` | VLP17, VLP20, VLP_100, bare gold NPs | Gold NP (near-black, circular) | `analysis/vlp_measure_v2.py` |
| `BMV` | BMV, BOG | Bright protein ring + dark stain pool (donut signature) | `analysis/bmv_measure.py` |

Sample classification from a raw image is **not yet automated** — the scientist supplies the sample type, or the agent infers it from filenames.

## CLI reference (all scripts, all flags)

Every script in `analysis/` is invocable directly. Canonical form is `uv run python -m analysis.<name>`.

### `analysis.vlp_measure_v2` — VLP gold + capsid measurement

```bash
uv run python -m analysis.vlp_measure_v2 <image_path> [flags]
```

| flag | default | what it does |
|---|---|---|
| `image_path` (positional) | required | A `.dm3`/`.dm4` file or a folder of them |
| `--pattern` | `*` | Glob when `image_path` is a folder |
| `--out` | `results/vlp_v2` | Output directory |
| `--sample-name` | folder name | Override the auto-derived sample name |
| `--workers` | `4` | Parallel worker processes |
| `--gold-threshold` | auto (Otsu) | Override the per-image gold-detection intensity cutoff |
| `--min-gold-nm` | `7.0` | Lower size cutoff for gold detection |
| `--max-gold-nm` | `30.0` | Upper size cutoff for gold detection |
| `--show-flagged` | off | Also draw flagged-unreliable circles in the overlay (tomato) |
| `--gui` | off | Pop a Tk progress window |

### `analysis.bmv_measure` — BMV / BOG capsid measurement

```bash
uv run python -m analysis.bmv_measure <image_path> [flags]
```

| flag | default | what it does |
|---|---|---|
| `image_path` (positional) | required | A `.dm3`/`.dm4` file or a folder of them |
| `--pattern` | `*` | Glob when `image_path` is a folder |
| `--out` | `results/bmv` | Output directory |
| `--expected-nm` | `28.0` | Expected capsid diameter (Hough searches around this) |
| `--workers` | `4` | Parallel worker processes |
| `--show-flagged` | off | Draw flagged (tomato) and no-fit (magenta) circles in overlay |
| `--debug` | off | Save radial profile plots for a sample of particles |
| `--debug-indices N1 N2 …` | random sample | Specific particle indices to plot in debug mode |
| `--gui` | off | Pop a Tk progress window |

### `analysis.eval` — compare a run to the reference benchmarks

```bash
uv run python -m analysis.eval <run_dir> [flags]
```

| flag | default | what it does |
|---|---|---|
| `run_dir` (positional) | required | Directory with the measurement CSV (e.g. `results/<sample>`) |
| `--sample-type` | `VLP` | `VLP` or `BMV` |
| `--no-report` | off | Skip writing `eval_report.md`; still returns the eval dict |
| `--json` | off | Print full result as JSON |

### `analysis.add_to_reference run` — append a run to `reference_runs.csv`

```bash
uv run python -m analysis.add_to_reference run <run_dir> --approver <name> [flags]
```

| flag | default | what it does |
|---|---|---|
| `run_dir` (positional) | required | Directory with the measurement CSV |
| `--sample-type` | `VLP` | `VLP` or `BMV` |
| `--sample-name` | folder name | Override the auto-derived sample name |
| `--approver` | required | Person approving this run (audit trail) |
| `--notes` | empty | Free-form note explaining why this run is reference-worthy |
| `--force` | off | Replace existing rows on `(sample_name, filename)` collision |

### `analysis.add_to_reference hand` — append hand measurements to `reference_hand.csv`

```bash
uv run python -m analysis.add_to_reference hand <hand_csv> --sample-name <name> --scientist <name> [flags]
```

| flag | default | what it does |
|---|---|---|
| `hand_csv` (positional) | required | ImageJ-exported CSV |
| `--sample-name` | required | Identifier shared with the script run this hand data validates |
| `--sample-type` | `VLP` | `VLP` or `BMV` |
| `--unit` | `um` | Length unit in source CSV: `um` or `nm` |
| `--scientist` | required | Person who hand-measured (audit trail) |
| `--measure-date` | today | ISO date when measurements were taken |
| `--run-dir` | none | Also copy the hand CSV into `<run-dir>/hand/` so per-batch eval picks it up |
| `--force` | off | Replace existing rows on `(sample_name, source_file)` collision |

### `analysis.seed_benchmarks` — initial seed (no flags)

```bash
uv run python -m analysis.seed_benchmarks
```

One-shot. Reads existing run directories under `results/` and writes `benchmarks/<sample_type>/{reference_runs,reference_hand}.csv`. Re-running overwrites — use `add_to_reference` for incremental growth.

### `analysis.plot_vlp_scatter` — combined scatter / 2D-hist / KDE / pooled hist

```bash
uv run python -m analysis.plot_vlp_scatter <csv> [<csv> …] [flags]
```

| flag | default | what it does |
|---|---|---|
| `csvs` (positional, ≥1) | required | One or more `vlp_measurements.csv` files |
| `--out-dir` | first CSV's folder | Where the combined plots land |
| `--combined-hist2d` | off | One pooled 2D histogram instead of per-sample subplots |

### `analysis.plot_capsid_groups` — per-image strip plot

```bash
uv run python -m analysis.plot_capsid_groups <csv> [flags]
```

| flag | default | what it does |
|---|---|---|
| `csv` (positional) | required | A `vlp_measurements.csv` |
| `--out` | CSV's folder | Where `capsid_groups.png` lands |

## Repo contents

| | |
|---|---|
| `analysis/` | Measurement, eval, plotting, reference-set management. |
| `benchmarks/<sample_type>/` | `reference_runs.csv` (script outputs we trust) + `reference_hand.csv` (per-particle hand measurements). |
| `incoming/` | DM3 dump location (gitignored in scientist repos). |
| `results/` | Per-run outputs (gitignored in scientist repos). |
| `LAB_NOTEBOOK.md` | Decision log: every non-trivial measurement decision, with dates and reasons. |
| `CLAUDE.md` | Co-scientist working principles + auto-bootstrap rule. Codex users: rename to `AGENTS.md`. |
| `.claude/skills/lab-pipeline/` | Project-scoped Claude Code skill. |
| `bootstrap.sh` | One-shot installer for new scientist repos. |

## Working principles

See `CLAUDE.md` for the full list. Non-negotiables:

1. **Always look at the actual image first** before picking a detector or diagnosing a bad result.
2. **Generate overlay PNGs for every run** — numbers alone cannot validate a pipeline.
3. **Auto-threshold per image.** CLAHE and staining variation mean fixed cutoffs misbehave.
4. **Declare physical scales in nm**, convert to pixels per-image.
5. **Surface anomalies** — 100% detection rates, step changes, outliers. Print and flag.

## Planned but not built

- `analysis/validate_script.py` — re-run the script on every image in `reference_runs.csv` after a script change; flag any image whose new measurement differs by more than tolerance from the recorded reference. The regression-test for the measurement pipeline itself.
