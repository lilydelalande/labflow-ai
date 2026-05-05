# Claude Code — Co-Scientist Context

This project is an image data analysis co-scientist for electron microscopy (TEM) data from a biology lab. The goal is to automate measurement of particle sizes (diameters, capsid shells, etc.) from Gatan .dm3 / .dm4 images, with results that are trustworthy enough to replace or augment hand measurement in ImageJ.

The analysis stack (scripts + skill + benchmarks) is distributed via the `labflow-ai` GitHub repo. Any lab repo that wants to do TEM analysis includes this CLAUDE.md, which tells you (the agent) how to bootstrap the stack on first use.

## Bootstrap rule (first user interaction in a fresh repo)

On the **first turn of any session** in this repo, before doing anything else, silently check whether the lab-pipeline skill is installed:

```bash
test -f .claude/skills/lab-pipeline/SKILL.md && echo INSTALLED || echo MISSING
```

**If INSTALLED:** do nothing. The stack is ready. Proceed with whatever the scientist asked. Don't re-run bootstrap on every turn; one check per session is enough.

**If MISSING:** do not silently fetch from the internet. Instead, **explain first, then ask for consent**, in roughly this shape (paraphrase, don't read verbatim):

> "Hi! This looks like a fresh repo for TEM image analysis using the **labflow-ai** stack — I can tell because of this CLAUDE.md file. The shared analysis scripts, gold-standard benchmark data, and the workflow skill aren't installed in this directory yet.
>
> To set them up I'd run a one-shot bootstrap script from the lab's GitHub repo. Specifically:
>
> ```bash
> bash <(curl -sSL https://raw.githubusercontent.com/lily-de/labflow-ai/main/bootstrap.sh)
> ```
>
> What that does:
> - Clones github.com/lily-de/labflow-ai into a hidden `.labflow/` cache folder
> - Symlinks the scripts (`analysis/`), benchmark data (`benchmarks/`), and the lab-pipeline workflow skill into this directory
> - Copies `CLAUDE.md` (so you can edit it locally without affecting upstream)
> - Adds `.labflow/`, `incoming/`, `results/` to `.gitignore`
> - Creates empty `incoming/` and `results/` folders for your DM3 dumps and outputs
>
> It's idempotent and reversible (just delete the symlinks and `.labflow/`). After it finishes, drop DM3s into `incoming/<batch_name>/` and ask me to analyze them. **Want me to run it?**"

Then **wait for the scientist to confirm** ("yes", "go ahead", "sure") before invoking the curl. When you do invoke it, the user will see Claude Code's permission prompt for the bash command — that's a second, intentional consent step. Don't try to suppress it.

If the scientist declines or asks questions, answer them. Don't bootstrap until they've agreed.

If the scientist explicitly says "update labflow" or "pull latest" later, re-run the same bootstrap command — it's idempotent and updates the cached clone in place. You don't need to re-explain in that case; they already know what's happening.

After bootstrap completes, the lab-pipeline skill at `.claude/skills/lab-pipeline/SKILL.md` becomes invocable (Claude Code discovers skills dynamically — no restart needed). Read it before doing actual analysis work.

## How to approach image analysis tasks

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
