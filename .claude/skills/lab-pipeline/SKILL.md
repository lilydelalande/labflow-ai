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
    image_path="incoming/2026-05-15_VLP17_batch4",   # file or directory
    sample_type="VLP",                               # "VLP" | "BMV" (BMV is a separate script today)
    batch_id=None,                                   # auto: folder name (or filename + today's date for single-image). Override only if scientist supplies a specific tag.
    out_dir="results/vlp17_batch4",
    pattern="VLP17_*",                               # glob without extension
    workers=6,                                       # parallelism
    show_flagged=False,                              # draw flagged circles in overlay
)
```

Returns:

```python
{
    "sample_type":    "VLP",
    "batch_id":       "VLP17_2026-05-15",
    "script_version": "vlp_measure_v2@2.0(wall=0.75;contrast=0.05;cv=0.2;...)",
    "summary": {
        "n_images": 17, "n_gold_total": 4060, "n_reliable": 2964,
        "n_dropped": 1096, "drop_rate": 0.27,
        "gold_mean_nm": 16.7, "gold_std_nm": 1.4, "gold_median_nm": 16.7,
        "capsid_mean_nm": 34.4, "capsid_std_nm": 1.8, "capsid_median_nm": 34.2,
    },
    "per_image": [
        {
            "filename": "VLP17_0001.dm4",
            "n_gold": 109, "n_wall_fit": 68, "n_reliable": 49,
            "wall_fit_success_rate": 0.62, "reliable_rate": 0.45, "drop_rate": 0.55,
            "gold_median_nm": 16.93, "gold_std_nm": 1.66,
            "capsid_median_nm": 32.94, "capsid_std_nm": 1.81,
            "median_wall_cv": 0.13, "iqr_wall_cv": 0.06,
        },
        ...
    ],
    "outputs": {
        "run_dir":         "results/vlp17_batch4",
        "csv_path":        "results/vlp17_batch4/vlp_measurements.csv",
        "overlays_dir":    "results/vlp17_batch4/overlays",
        "histograms_path": "results/vlp17_batch4/vlp_histograms.png",
        "scatter_path":    "results/vlp17_batch4/vlp_scatter.png",
        "summary_md":      "results/vlp17_batch4/SUMMARY.md",
    },
}
```

The same data is in `SUMMARY.md` as a human-readable markdown table. Read it to the user, link the overlays and plots — don't re-derive numbers from the CSV when the dict already has them.

### Calling pattern

Prefer `Bash`:

```bash
uv run python -c "
import json
from analysis.vlp_measure_v2 import run
r = run(image_path='...', batch_id='...', out_dir='...', workers=6)
print(json.dumps(r, indent=2))
"
```

The full result is small (kilobytes) — fine to ingest as one tool result. If a user wants to run the same thing via the terminal, they invoke `uv run python -m analysis.vlp_measure_v2 …` directly; the CLI wrapper calls the same `run()` and prints a headline.

`run()` automatically calls `evaluate()` at the end (if `benchmarks/vlp/reference_runs.csv` exists), so a single call gives you both `SUMMARY.md` and `eval_report.md`. Inspect `result["eval"]` for the headline; full report is at `result["outputs"]["eval_report"]`.

### `eval.evaluate(run_or_dir, sample_type="VLP")` — compare against the reference

Two checks per the LAB_NOTEBOOK design:

- **Per-image quality check** — flags any image whose dimensionless metrics (`wall_fit_success_rate`, `reliable_rate`, `drop_rate`, `median_wall_cv`, `capsid_*`) sit outside median ± 2 × IQR of the same-subtype reference distribution.
- **Hand vs script** — if `reference_hand.csv` has rows tagged with this batch_id (or, fallback, the same sample_subtype), reports `hand_capsid_mean − script_capsid_mean`.

Accepts either a `run()` dict or a path to a run directory. Auto-called inside `run()`; rarely needed standalone unless re-evaluating an old run after the reference grows.

### `add_to_reference.add_run(...)` / `add_to_reference.add_hand(...)` — grow the reference

The **manual gate** for promoting new data into `benchmarks/vlp/`. Both functions:

- **Refuse duplicates by default.** `add_run` keys on `(batch_id, filename)`; `add_hand` keys on `(batch_id, source_file)`. If the data is already in the reference, you get a `DuplicateReferenceError` listing the colliding rows. Pass `force=True` to *replace* the existing rows (e.g. when re-running a script version).
- **Require `approver` (run) or `scientist` (hand)** as a non-empty string. Recorded in the CSV as the audit trail.

```python
from analysis.add_to_reference import add_run, add_hand

add_run(
    run_dir       = "results/2026-05-15_VLP17_batch4",
    sample_type   = "VLP",
    approver      = "Lily",
    notes         = "Hand-validated; replaces prior VLP17 reference.",
    # batch_id auto-derives from run_dir folder name; sample_subtype auto-infers per image
)

add_hand(
    hand_csv     = "results/2026-05-15_VLP17_batch4/hand/measurements.csv",
    batch_id     = "2026-05-15_VLP17_batch4",   # match the run's batch_id so eval can pair them
    sample_type  = "VLP",
    hand_format  = "paired",          # or "capsid_only" for per-image diagnostic CSVs
    length_unit  = "um",              # or "nm" — depends on ImageJ calibration
    scientist    = "Lily",
    measure_date = "2026-05-15",
)
```

When to invoke `add_run` / `add_hand`:
- The scientist explicitly asks: "add this run to the reference", "save these as benchmarks", "this batch should be the new baseline", etc.
- Never invoke automatically. Promotion to reference is always a deliberate human decision.

When the dedup error fires, surface the exact collision to the scientist (the message lists the colliding `batch_id` + `filename`). Don't auto-`force` — let the scientist decide whether the existing rows should be replaced.

## File conventions

All paths are relative to the **current working directory** (the repo root the scientist is in). Do not write outside it.

- `incoming/<batch_name>/` — DM3/DM4 dump location, one subfolder per batch. Filenames inside the folder are arbitrary.
- `results/<batch_name>/` — measurement output: CSV, overlays, plots, SUMMARY.md. Mirrors the incoming folder name 1:1.
- `results/<batch_name>/hand/` — hand-measurement CSVs for this batch (when available).
- `benchmarks/<sample_type>/gold_standard.csv` — curated approved runs (manual gate). Not yet implemented; design banked in `LAB_NOTEBOOK.md`.
- `benchmarks/<sample_type>/hand_measurements.csv` — per-particle hand data, joined to runs by `batch_id`. Not yet implemented.
- `LAB_NOTEBOOK.md` — append a dated section any time you make a non-trivial decision (new sample type, threshold change, etc).

### First-time directory setup
If `incoming/`, `results/`, or `benchmarks/` don't exist in the current working directory, **create them on first use** — but explain what you're doing first, in one short message:

> "I don't see `incoming/` or `results/` here yet. I'll create them in the current directory. New batches go in `incoming/<batch_name>/`, results land in `results/<batch_name>/`, and gold-standard reference data lives in `benchmarks/<sample_type>/`."

Then run `mkdir -p incoming results benchmarks` and continue. Don't ask permission — these are inert empty folders, fully reversible. Do this once per fresh repo, never again.

## Standard scientist workflows

### "Analyze the new images" / "measure batch X"

The scientist should never have to specify `sample_type` if you can figure it out yourself.

1. **Locate the batch folder.** If they named it (`"the VLP17 batch I just dropped"`), look in `incoming/` for the obvious match. If they didn't, list `incoming/*` and ask only if it's ambiguous.
2. **Infer `sample_type` from the folder/filename.** Folder names like `VLP17_*`, `VLP_100_*`, `VLP20_*` → `sample_type="VLP"`. Names containing `BMV` or `BOG` → `sample_type="BMV"`. If the folder name is uninformative, peek at one of the filenames inside (`VLP17_0001.dm4` → VLP). Only ask the scientist if both folder name and filenames are uninformative.
3. **`batch_id` auto-derives** from the folder name — don't pass it explicitly.
4. **Set `out_dir=results/<folder_name>`** so results mirror the incoming folder structure.
5. **Call `run(...)`** and read the resulting `SUMMARY.md`.
6. **Quote the headline back** — n_reliable, capsid mean ± std, drop rate, link to overlays + SUMMARY.md.
7. **Flag suspicious images:** any `per_image` entry with `reliable_rate < 0.7` or `wall_fit_success_rate < 0.85` is a candidate for hand validation. Name them explicitly in the response.

The scientist's bar should be: drop a folder in `incoming/`, say "analyze it" with at most a folder hint. The agent does the rest. Asking "what sample type?" when it's `VLP17_*.dm4` files is a failure mode — don't.

### "Compare to my hand measurements"
1. Add the hand CSV to the reference (with the same `batch_id` as the run): `add_hand(...)`.
2. Re-evaluate the run: `evaluate(run_dir)` — Hand vs script section will now populate.
3. Quote the resulting delta back. If it sits well outside historical hand-vs-script deltas in `reference_hand.csv` for that subtype, flag it.

### "Add this run to the reference"
1. Confirm the scientist actually wants this run promoted (it's a deliberate decision, not a default).
2. Call `add_to_reference.add_run(run_dir, sample_type, approver, notes)`.
3. If `DuplicateReferenceError` fires, surface the colliding rows verbatim. Ask the scientist whether to replace (`force=True`) or skip.

### "Did the script regress?"
Re-run measurement on the images already represented in `reference_runs.csv` (group by `batch_id`, find the source images, re-run). Then compare new `per_image` rows to the reference rows for those images. Differences > ~0.3 nm in capsid_median are regressions worth surfacing.

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
