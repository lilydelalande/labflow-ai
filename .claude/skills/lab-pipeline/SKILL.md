---
name: lab-pipeline
description: Drive the TEM image-analysis pipeline — measure VLP/BMV samples, eval against the gold-standard set, integrate hand measurements. Use this skill when the scientist asks to analyze, measure, or evaluate TEM (Gatan .dm3/.dm4) data.
---

# Lab pipeline — TEM measurement and eval

This skill is the canonical interface for running the TEM image-analysis pipeline. The codebase exposes pure Python functions (no shell scraping) so you should **import and call them directly**, then quote the structured return values back to the user. Avoid running the CLI wrappers via Bash unless the user explicitly asks for a terminal-style run.

## When to invoke this skill

The scientist asks any of:
- "Measure this batch of DM3s" / "analyze this folder"
- "Run the VLP/BMV measurement script on…"
- "How does this batch compare to our gold-standard?"
- "Compare the script output to my hand measurements"
- "Add this run to the gold standard"
- "Did the script regress on the gold-standard set?"

If the request is exploratory ("show me the radial profile of particle 12") use the existing scripts directly — this skill is for the standard measure → eval → report workflow.

## Tool surface

Three importable Python tools, all in `analysis/`:

1. `vlp_measure_v2.run(...)` — measurement (always; auto-runs eval if reference exists)
2. `eval.evaluate(...)` — compare a run to the reference benchmarks
3. `add_to_reference.add_run(...)` / `add_to_reference.add_hand(...)` — manually grow the reference set

Prefer importing and calling these directly over the CLI wrappers — the structured returns are the cleanest interface.

### `vlp_measure_v2.run()` — measure a VLP sample

Located in `analysis/vlp_measure_v2.py`. Call from Python:

```python
from analysis.vlp_measure_v2 import run

result = run(
    image_path="incoming/<sample_name>",   # file or directory
    sample_type="VLP",                     # "VLP" | "BMV"
    sample_name=None,                      # auto: input folder name. Override only if scientist asks.
    out_dir="results/<sample_name>",       # mirrors the input folder by convention
    pattern="*",                           # glob without extension
    workers=6,
    show_flagged=False,
)
```

Returns:

```python
{
    "sample_type":    "VLP",
    "sample_name":    "<input folder name>",
    "script_version": "vlp_measure_v2@2.0(...)",
    "summary":        { n_images, n_gold_total, n_reliable, gold_*, capsid_*, ... },
    "per_image":      [ {filename, n_gold, n_reliable, wall_fit_success_rate, ..., capsid_median_nm, median_wall_cv, iqr_wall_cv}, ... ],
    "outputs":        { run_dir, csv_path, overlays_dir, histograms_path, scatter_path, summary_md, eval_report (if benchmarks exist) },
    "eval":           { n_warns_total, n_ref_runs, hand_vs_script, report_path },  # auto-runs after measurement
}
```

The same headline is in `SUMMARY.md` (human-readable). Read it to the user, link the overlays — don't re-derive numbers from the CSV when the dict already has them.

### Calling pattern

```bash
uv run python -c "
import json
from analysis.vlp_measure_v2 import run
r = run(image_path='incoming/<sample_name>', out_dir='results/<sample_name>', workers=6)
print(json.dumps(r, indent=2))
"
```

CLI alternative for terminal use: `uv run python -m analysis.vlp_measure_v2 …`.

`run()` automatically calls `evaluate()` at the end (if `benchmarks/vlp/reference_runs.csv` exists), so a single call gives you both `SUMMARY.md` and `eval_report.md`.

### `eval.evaluate(run_or_dir, sample_type="VLP")` — compare against the reference

Two checks:

- **Per-image quality check** — for each image, compares its **size-invariant quality metrics** (`wall_fit_success_rate`, `reliable_rate`, `drop_rate`, `median_wall_cv`, `iqr_wall_cv`) against the pooled distribution of the same metrics across all prior trusted runs in `reference_runs.csv`. Flags values outside median ± 2 × IQR. Absolute sizes (`capsid_mean_nm`, `gold_mean_nm`) are reported but not gated — sample names are intent, not measured size, and the eval shouldn't lock to mislabelled bins.
- **Hand vs script** — if any `*.csv` exist in `<run_dir>/hand/`, parse them (paired gold/capsid format, ImageJ Length column), compute hand_mean − script_mean for capsid and gold. Surfaces calibration drift between the script and human ImageJ measurements. Filesystem-paired only — no batch IDs, no fuzzy subtype matching.

Auto-called inside `run()`; standalone use is for re-evaluating an old run after the reference grows.

### `add_to_reference.add_run(...)` / `add_to_reference.add_hand(...)` — grow the reference

The **manual gate** for promoting new data into `benchmarks/vlp/`. Both functions:

- **Refuse duplicates by default.** `add_run` dedups on `(sample_name, filename)`; `add_hand` dedups on `(sample_name, source_file)`. Pass `force=True` to replace existing rows.
- **Require `approver` (run) or `scientist` (hand)** as a non-empty string — audit trail recorded in the CSV.

```python
from analysis.add_to_reference import add_run, add_hand

add_run(
    run_dir     = "results/<sample_name>",
    sample_type = "VLP",
    approver    = "Lily",
    notes       = "hand-validated, looks clean",
    # sample_name auto-derives from run_dir folder name
)

add_hand(
    hand_csv    = "path/to/imagej_export.csv",
    sample_name = "<sample_name>",        # matches the run's sample_name so eval can pair them
    sample_type = "VLP",
    length_unit = "um",                   # or "nm" — match what ImageJ exported
    scientist   = "Lily",
    run_dir     = "results/<sample_name>", # also copies CSV into <run_dir>/hand/ for filesystem-paired eval
)
```

When to invoke `add_run` / `add_hand`:
- The scientist explicitly asks: "add this run to the reference", "save these as benchmarks", "this batch should be the new baseline", etc.
- Never invoke automatically. Promotion to reference is always a deliberate human decision.

When the dedup error fires, surface the exact collision to the scientist (the error lists the colliding `(sample_name, filename)` or `(sample_name, source_file)` tuples). Don't auto-`force` — let the scientist decide.

## File conventions

All paths are relative to the **current working directory** (the repo root the scientist is in). Do not write outside it.

- `incoming/<sample_name>/` — DM3/DM4 dump location. The folder name IS the `sample_name` identifier — scientist names it whatever they want.
- `results/<sample_name>/` — measurement output: CSV, overlays, plots, SUMMARY.md, eval_report.md. Mirrors the incoming folder name 1:1.
- `results/<sample_name>/hand/` — hand-measurement CSVs for this run (when available). Filesystem-paired to the run; eval reads any `*.csv` in here.
- `benchmarks/<sample_type>/reference_runs.csv` — curated approved runs (manual gate via `add_to_reference.add_run`).
- `benchmarks/<sample_type>/reference_hand.csv` — accumulated hand measurements (manual gate via `add_to_reference.add_hand`).
- `LAB_NOTEBOOK.md` — append a dated section any time you make a non-trivial decision (new sample type, threshold change, etc).

### First-time directory setup
If `incoming/`, `results/`, or `benchmarks/` don't exist in the current working directory, **create them on first use** — but explain what you're doing first, in one short message:

> "I don't see `incoming/` or `results/` here yet. I'll create them in the current directory. New samples go in `incoming/<sample_name>/`, results land in `results/<sample_name>/`, and reference data lives in `benchmarks/<sample_type>/`."

Then run `mkdir -p incoming results benchmarks` and continue. Don't ask permission — these are inert empty folders, fully reversible. Do this once per fresh repo, never again.

## Standard scientist workflows

### "Analyze the new images" / "measure batch X"

The scientist should never have to specify `sample_type` if you can figure it out yourself.

1. **Locate the sample folder.** If they named it (`"the VLP17 sample I just dropped"`), look in `incoming/` for the obvious match. If they didn't, list `incoming/*` and ask only if it's ambiguous.
2. **Infer `sample_type` from the folder/filename.** Folder names like `VLP17_*`, `VLP_100_*`, `VLP20_*` → `sample_type="VLP"`. Names containing `BMV` or `BOG` → `sample_type="BMV"`. If the folder name is uninformative, peek at one of the filenames inside (`VLP17_0001.dm4` → VLP). Only ask the scientist if both folder name and filenames are uninformative.
3. **`sample_name` auto-derives** from the input folder name — don't pass it explicitly.
4. **Set `out_dir=results/<folder_name>`** so results mirror the incoming folder structure.
5. **Call `run(...)`** and read the resulting `SUMMARY.md`.
6. **Quote the headline back** — n_reliable, capsid mean ± std, drop rate, link to overlays + SUMMARY.md.
7. **Flag suspicious images:** any `per_image` entry with `reliable_rate < 0.7` or `wall_fit_success_rate < 0.85` is a candidate for hand validation. Name them explicitly in the response.

The scientist's bar should be: drop a folder in `incoming/`, say "analyze it" with at most a folder hint. The agent does the rest. Asking "what sample type?" when it's `VLP17_*.dm4` files is a failure mode — don't.

### "Compare to my hand measurements"
1. The simplest path: drop the ImageJ CSV directly into `results/<sample_name>/hand/`. Re-evaluate (`evaluate(run_dir)`); the Hand vs script section populates automatically.
2. To also save the hand CSV into the reference for future runs, call `add_to_reference.add_hand(...)` with `run_dir=results/<sample_name>` — that copies it into `hand/` AND appends rows to `reference_hand.csv`.
3. Quote the delta back to the scientist. A small delta (a few tenths of a nm) is normal noise. A persistent positive or negative delta across multiple runs would suggest the script is calibrated differently from the human's tracing — worth investigating.

### "Add this run to the reference"
1. Confirm the scientist actually wants this run promoted (it's a deliberate decision, not a default).
2. Call `add_to_reference.add_run(run_dir, sample_type, approver, notes)`.
3. If `DuplicateReferenceError` fires, surface the colliding rows verbatim. Ask the scientist whether to replace (`force=True`) or skip.

### "Did the script regress?"
Re-run measurement on the images already represented in `reference_runs.csv` (group by `sample_name`, find the source images, re-run). Then compare new `per_image` rows to the reference rows for those images. Differences > ~0.3 nm in capsid_median are regressions worth surfacing.

## Output discipline

Three rules — keep them tight:

1. **Every tool writes both structured JSON (return value) AND human artifacts to disk.** The agent reads the dict; the human reads the markdown / plots / overlays. Don't choose one.
2. **Quote the headline back.** When you finish a `run()` call, surface: n_reliable, capsid mean ± std, drop rate, link to SUMMARY.md and overlays. Don't make the scientist hunt for it in the conversation.
3. **Never silently swallow flagged images.** If `per_image` shows entries with very low `reliable_rate`, name them in the response — that's the kind of thing the scientist needs to see.

## Sample-type routing

| sample_type | script | morphology |
|---|---|---|
| `VLP` (any gold-anchored sample: VLP17, VLP20, VLP_100) | `vlp_measure_v2.run()` | gold core + protein wall |
| `BMV`, `BOG` (no gold) | `bmv_measure.py` (CLI only for now) | dark-stain center + bright protein ring |

Sample classification from a raw image is **not yet automated** — always ask the scientist for the sample type if `sample.txt` doesn't specify.

## What this skill does NOT do

- Auto-classify a DM3 — ask the scientist.
- Decide whether a run is "good enough" to be promoted to gold — that is an explicit human gate.
- Touch `benchmarks/` files automatically. All gold-standard mutations are explicit and recorded in git.
- Run the BMV script via this skill yet — that's planned but not Phase-1.
