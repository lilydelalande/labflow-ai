#!/usr/bin/env bash
# bootstrap.sh — install the labflow-ai analysis stack into the current directory.
#
# Idempotent: safe to re-run. Updates the cached clone in place; rebuilds
# symlinks; preserves the local CLAUDE.md (which scientists may edit).
#
# Usage:
#   curl -sSL https://raw.githubusercontent.com/lily-de/labflow-ai/main/bootstrap.sh | sh
#   bash <(curl -sSL https://raw.githubusercontent.com/lily-de/labflow-ai/main/bootstrap.sh)
#   ./bootstrap.sh --relink     # only refresh symlinks; don't pull
#
# What it does:
#   1. Clones (or pulls) labflow-ai into .labflow/
#   2. Symlinks:
#        analysis/                  -> .labflow/analysis
#        benchmarks/                -> .labflow/benchmarks
#        .claude/skills/lab-pipeline -> .labflow/.claude/skills/lab-pipeline
#   3. Copies CLAUDE.md if not already present
#   4. Adds .labflow/, incoming/, results/ to .gitignore
#   5. Records the upstream SHA in .labflow/INSTALLED_SHA
#   6. Creates incoming/ and results/ scaffolding folders

set -euo pipefail

REPO_URL="${LABFLOW_REPO:-https://github.com/lily-de/labflow-ai.git}"
BRANCH="${LABFLOW_BRANCH:-main}"
CACHE_DIR=".labflow"

mode="install"
if [ "${1:-}" = "--relink" ]; then
    mode="relink"
fi

say() { printf "  %s\n" "$*"; }
hdr() { printf "\n→ %s\n" "$*"; }

hdr "labflow-ai bootstrap"
say "repo: $REPO_URL ($BRANCH)"
say "cwd:  $(pwd)"

# 1. Clone or update the cache
if [ "$mode" = "install" ]; then
    if [ -d "$CACHE_DIR/.git" ]; then
        hdr "Updating existing $CACHE_DIR"
        git -C "$CACHE_DIR" fetch --depth 1 origin "$BRANCH"
        git -C "$CACHE_DIR" reset --hard "origin/$BRANCH"
    else
        hdr "Cloning $REPO_URL into $CACHE_DIR"
        git clone --depth 1 --branch "$BRANCH" "$REPO_URL" "$CACHE_DIR"
    fi
fi

if [ ! -d "$CACHE_DIR" ]; then
    echo "ERROR: $CACHE_DIR does not exist. Run without --relink first." >&2
    exit 1
fi

INSTALLED_SHA="$(git -C "$CACHE_DIR" rev-parse HEAD)"
echo "$INSTALLED_SHA" > "$CACHE_DIR/INSTALLED_SHA"
say "pinned to $INSTALLED_SHA"

# 2. Symlinks
hdr "Linking scripts, benchmarks, and skill"
ln -sfn "$CACHE_DIR/analysis"   analysis
say "analysis/ -> $CACHE_DIR/analysis"

ln -sfn "$CACHE_DIR/benchmarks" benchmarks
say "benchmarks/ -> $CACHE_DIR/benchmarks"

mkdir -p .claude/skills
ln -sfn "../../$CACHE_DIR/.claude/skills/lab-pipeline" .claude/skills/lab-pipeline
say ".claude/skills/lab-pipeline -> $CACHE_DIR/.claude/skills/lab-pipeline"

# 3. Copy CLAUDE.md (only if absent — preserve local edits)
if [ ! -f CLAUDE.md ]; then
    cp "$CACHE_DIR/CLAUDE.md" CLAUDE.md
    say "CLAUDE.md copied (you can edit it locally; commit upstream PRs to share changes)"
else
    say "CLAUDE.md already exists — leaving it alone"
fi

# 4. .gitignore
hdr "Updating .gitignore"
touch .gitignore
for entry in ".labflow/" "incoming/" "results/"; do
    if ! grep -qxF "$entry" .gitignore; then
        printf "%s\n" "$entry" >> .gitignore
        say "added $entry"
    else
        say "$entry already ignored"
    fi
done

# 5. Scaffold working directories
hdr "Scaffolding incoming/ and results/"
mkdir -p incoming results

# 6. Done
hdr "Ready"
cat <<EOF
  Drop DM3/DM4 files into incoming/<batch_name>/ and either:

    Talk to Claude / Codex:  "analyze the new images in incoming/<batch_name>"
    Run scripts directly:    uv run python -m analysis.vlp_measure_v2 incoming/<batch_name> --workers 6

  Both paths write to results/<batch_name>/ — identical outputs.

  Update later:  bash <(curl -sSL https://raw.githubusercontent.com/lily-de/labflow-ai/main/bootstrap.sh)

EOF
