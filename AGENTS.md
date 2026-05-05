# Codex — Co-Scientist Context

This project is an image data analysis co-scientist for electron microscopy (TEM) data from a biology lab. The goal is to automate measurement of particle sizes (diameters, capsid shells, etc.) from Gatan .dm3 / .dm4 images, with results that are trustworthy enough to replace or augment hand measurement in ImageJ.

The analysis stack (scripts + skill + benchmarks) is distributed via the `labflow-ai` GitHub repo. Any lab repo that wants to do TEM analysis includes this AGENTS.md (and CLAUDE.md), which tells you (the agent) how to bootstrap the stack on first use.

## Bootstrap rule (run this check on every first user interaction)

On the **first turn of any session** in this repo, before doing anything else, check whether the lab-pipeline skill is installed:

```bash
test -f .claude/skills/lab-pipeline/SKILL.md && echo INSTALLED || echo MISSING
```

If MISSING, immediately tell the scientist (in one short sentence) that you're bootstrapping, then run:

```bash
bash <(curl -sSL https://raw.githubusercontent.com/lily-de/labflow-ai/main/bootstrap.sh)
```

This is idempotent and safe. It clones labflow-ai into `.labflow/`, symlinks `analysis/`, `benchmarks/`, and the lab-pipeline skill into the working directory, and gitignores `.labflow/`, `incoming/`, `results/`. Wait for it to finish, then proceed with whatever the scientist actually asked.

If INSTALLED, do nothing — the stack is ready. Don't re-run bootstrap on every turn; one check at session start is enough.

The scientist should never have to think about installation. Bootstrap is a side effect of starting work, not something they ask for.

If the scientist explicitly says "update labflow" or "pull latest", run the bootstrap command again — it's idempotent and updates the cached clone in place.

## Workflow guide

After bootstrap, **read `.claude/skills/lab-pipeline/SKILL.md`** before doing any TEM analysis. That file is the canonical workflow guide — it documents the tool surface (`run()` from `analysis.vlp_measure_v2`), file conventions, sample-type routing, output discipline, and the standard "scientist asks X → do Y" patterns.

The path is `.claude/skills/...` because that's where Claude Code looks for skills automatically. The content is plain markdown with YAML frontmatter and is agent-neutral — Codex should just `cat` or `Read` it like any other doc. Don't be confused by the YAML header (`name:`, `description:` lines): that's only used by Claude Code's slash-command discovery; the body is what matters.

Treat any guidance in SKILL.md as authoritative for "how to drive the pipeline." This AGENTS.md is for the higher-level co-scientist principles (anchor on easy features, generate overlays, surface anomalies); the skill file is for the concrete how-to.

## How to approach image analysis tasks

**Decompose — anchor on the easiest feature first.**
When a sample has multiple measurable features, detect the highest-contrast, most unambiguous feature first and use it as an anchor for harder measurements. Never try to detect everything in one step. Always ask: what is the easiest thing to reliably detect in this image, and what can I build on top of it?

**Generate overlay images automatically for every run.**
Numbers alone cannot validate a pipeline. Always produce an overlay PNG showing detected boundaries on the raw image. The scientist uses this to immediately judge whether detections are correct. This is non-negotiable — never output only a CSV.

**Link diagnostics to each other explicitly.**
If you save radial profile plots, number the particles in the overlay to match. If you highlight a subset in a debug view, make those same particles visually distinct in the overlay. The scientist must be able to say "particle #11 looks wrong — show me its profile" and find it immediately.

**The scientist's visual judgment is ground truth during development.**
When the scientist says "this looks wrong," believe them and investigate. When they identify specific particles that look correct, use those as reference cases to understand what the algorithm should be doing. Your job is to make visual inspection as efficient as possible, not to replace it.

**Auto-threshold per image. Never use a fixed threshold across images.**
CLAHE and staining variation mean the same intensity cutoff will behave differently on different images. Use Otsu or equivalent per-image automatic threshold selection whenever possible. When you add a tunable parameter, ask: can this be derived automatically from the image?

**Generate multiple views of output data.**
Histograms, scatter plots, and per-image summary tables reveal different artifacts. Quantization shows up in scatter plots but not histograms. Per-image drift shows up in tables but not aggregate summaries. Always produce all three.

**Flag suspicious aggregate metrics — don't silently accept them.**
A 100% detection rate, a step change between sequential images, or an outlier in the range should be printed and flagged. It may be real biology or it may be an artifact. Surface it and let the scientist decide.

**Sub-pixel refinement matters for scatter plots.**
Integer pixel measurements produce quantization banding in scatter plots. Use parabolic interpolation around detected minima/edges to get continuous values.

**Use `uv` for all package management.**
Run scripts with `uv run python`. Add dependencies with `uv add`.

## Sample types and scripts

| Sample | Script | Anchors on |
|--------|--------|------------|
| VLPs with gold NP core | `analysis/vlp_measure_v2.py` | Gold NP (near-black, circular) |
| Bare gold NPs | `analysis/vlp_measure_v2.py` (gold only) | Gold NP |
| Plain viruses (e.g. BMV, BOG) | `analysis/bmv_measure.py` | Bright protein ring + dark stain pool |

## Evaluation approach

For any new sample type or script: measure a subsample by hand in ImageJ first, run the script, compare mean ± std and detection rate. Build trust before running blind. This is the lab's validation protocol.
