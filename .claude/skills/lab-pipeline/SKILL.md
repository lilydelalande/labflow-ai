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

## Tool surface (Phase 1 — measurement only)

The eval / approve / validate tools are still being designed. Right now the working tool is:

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
Phase 1 doesn't have a paired-comparison tool yet. Until it does:
1. Run `run()` if not already done.
2. Read the hand CSV. Compute per-image (or per-batch) hand mean ± std and compare to script `capsid_mean_nm`.
3. Report `delta = hand_mean − script_mean` and how it stacks up to prior comparisons in `LAB_NOTEBOOK.md`.

### "Add this run to the gold standard"
Not yet implemented. Direct the user to design choices documented in the lab notebook ("Implications for the data-analysis agent" section) and propose appending a row by hand to `benchmarks/<sample_type>/gold_standard.csv` when that file exists.

### "Did the script regress?"
Not yet implemented. Until then, eyeball: re-run on a representative subset of `benchmarks/<sample_type>/` and confirm the per-image numbers match the gold-standard CSV within ~0.2 nm.

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
