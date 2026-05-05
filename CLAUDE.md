# Claude Code — Co-Scientist Context

This project is an image data analysis co-scientist for electron microscopy (TEM) data from a biology lab. The goal is to automate measurement of particle sizes (diameters, capsid shells, etc.) from Gatan .dm3 / .dm4 images, with results that are trustworthy enough to replace or augment hand measurement in ImageJ.

The analysis stack (scripts + skill + benchmarks) is distributed via the `labflow-ai` GitHub repo. Any lab repo that wants to do TEM analysis includes this CLAUDE.md, which tells you (the agent) how to bootstrap the stack on first use.

## Bootstrap rule (first user interaction in a fresh repo)

On the **first turn of any session** in this repo, before doing anything else, follow this sequence:

**1. Greet and announce the check** (one sentence, conversational, no jargon):

> "Hi! Quick check first — making sure the lab's analysis tools are set up in this folder."

**2. Run the check:**

```bash
test -f .claude/skills/lab-pipeline/SKILL.md && echo INSTALLED || echo MISSING
```

**3a. If INSTALLED:** confirm briefly and ask what they want to do:

> "All set. What would you like to analyze?"

Then proceed normally. Don't re-run the check on subsequent turns; one per session is enough.

**3b. If MISSING:** ask permission to install, briefly:

> "Looks like the tools aren't installed yet — takes about 30 seconds the first time. Want me to set them up?"

If they say yes, run the bootstrap and proceed. If they ask what it does, answer plainly: it pulls the analysis scripts and benchmark data from the lab's GitHub repo, installs the Python packages the scripts need, and sets up `incoming/` and `results/` folders. Most of the time is the package install (~30s first run; near-instant on re-runs). It can be removed by deleting `.labflow/` and `.venv/`. Avoid terms like "idempotent", "symlinks", "clone" unless the scientist is technical and asks. If they decline or ask other questions, answer them and don't bootstrap until they've agreed.

The bootstrap command:

```bash
bash <(curl -sSL https://raw.githubusercontent.com/lily-de/labflow-ai/main/bootstrap.sh)
```

When you invoke it, the user sees Claude Code's permission prompt — that's a second consent step. Don't try to suppress it.

If the scientist explicitly says "update labflow" or "pull latest" later, re-run the same bootstrap command — it's idempotent and updates the cached clone in place. You don't need to re-explain in that case; they already know what's happening.

After bootstrap completes, the lab-pipeline skill at `.claude/skills/lab-pipeline/SKILL.md` becomes invocable (Claude Code discovers skills dynamically — no restart needed). Read it before doing actual analysis work.

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
