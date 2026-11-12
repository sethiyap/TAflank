#!/bin/bash
#PBS -P rg42
#PBS -q copyq
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/rg42+scratch/rg42
#PBS -l wd
#PBS -N viz_setup
#PBS -j oe
#
# One-time setup for PANGRAPH_ANNOTATED/PANGRAPH_CONTEXT_SUMMARY on Gadi:
# pull a container with matplotlib+scipy. Same reasoning and the same
# pattern as env/gadi_bakta_setup.sh -- no module on Gadi has been confirmed
# to carry matplotlib, and compute nodes have no outbound network to pull
# one at run time, so it has to be fetched once, here, on copyq.
#
#   cd /g/data/rg42/TAflank
#   qsub env/gadi_viz_setup.sh             # via copyq (recommended)
#   bash env/gadi_viz_setup.sh             # or directly on a login node
#
# Idempotent: re-running skips whatever is already pulled.
#
# WHICH IMAGE. jupyter/scipy-notebook is a well-known, officially maintained
# public image that bundles matplotlib+scipy+numpy+pandas with python3 on
# PATH -- picked because it needs nothing beyond `singularity pull`, the
# same zero-build approach env/gadi_bakta_setup.sh already uses, rather than
# a from-scratch Dockerfile this environment has no way to build or host.
# It is NOT a minimal image (multiple GB) -- fine on Gadi's storage, same
# tradeoff as Bakta's own multi-GB database.
#
# This is a first reasonable choice, not a verified one: nobody has run this
# script yet. Step 3 below smoke-tests the actual imports this project
# needs (matplotlib, scipy) before declaring success, and step 4 smoke-tests
# pangraph_blocks.py itself against a real graph.json if one is available.
# If either fails, override with VIZ_IMAGE=docker://some/other:tag and
# rerun -- the script is idempotent and will just re-pull.
set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
VIZROOT="${VIZROOT:-/g/data/rg42/viz_env}"
VIZ_IMAGE="${VIZ_IMAGE:-docker://jupyter/scipy-notebook:latest}"
TAG=$(basename "$VIZ_IMAGE" | tr ':' '-')
SIF="$VIZROOT/${TAG}.sif"

module purge 2>/dev/null || true
module load singularity
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-/scratch/rg42/$USER/singularity_cache}"
mkdir -p "$VIZROOT" "$SINGULARITY_CACHEDIR"

echo "image         : $VIZ_IMAGE"
echo "install into  : $SIF"
echo "singularity   : $(singularity --version)"
echo

echo "=== 1. container ==="
if [ -s "$SIF" ]; then
    echo "  already pulled: $SIF"
else
    singularity pull "$SIF" "$VIZ_IMAGE"
    echo "  pulled: $SIF  ($(du -h "$SIF" | cut -f1))"
fi
VIZ=(singularity exec -B /g/data/rg42,/scratch/rg42 "$SIF")
echo

echo "=== 2. python ==="
"${VIZ[@]}" python3 --version
echo

echo "=== 3. import check ==="
"${VIZ[@]}" python3 -c "
import matplotlib, scipy
print('matplotlib', matplotlib.__version__)
print('scipy     ', scipy.__version__)
" || { echo "  FAILED: this image does not have matplotlib+scipy importable" >&2
       echo "  try a different image: VIZ_IMAGE=docker://... $0" >&2; exit 1; }
# openpyxl is needed only by scripts/pangraph_context_summary.py (the Excel
# output), not by PANGRAPH_BLOCKS/PANGRAPH_ANNOTATED -- checked separately
# and WARNED rather than failing the whole setup over it, since it is not
# required for the core pipeline to run.
"${VIZ[@]}" python3 -c "import openpyxl; print('openpyxl  ', openpyxl.__version__)" \
    || echo "  WARNING: openpyxl not importable in this image -- scripts/pangraph_context_summary.py's Excel output will fail. Everything else (figures, tables) is unaffected."
echo

echo "=== 4. smoke test on one real graph, if available ==="
GRAPH_JSON=$(ls "$REPO"/results/pangraph/graph.json \
                "$REPO"/run_all*/results/pangraph/graph.json 2>/dev/null | head -1 || true)
if [ -z "$GRAPH_JSON" ]; then
    echo "  skipped: no results/pangraph/graph.json found yet (run the pipeline first)"
else
    T="${PBS_JOBFS:-/tmp}/viz_smoke"
    rm -rf "$T"; mkdir -p "$T"
    "${VIZ[@]}" python3 "$REPO/scripts/pangraph_blocks.py" \
        --graph "$GRAPH_JSON" --outdir "$T" --out-tables "$T" \
        > "$T/log" 2>&1 \
        && echo "  OK  $(ls "$T"/pangraph_blocks.png 2>/dev/null && echo 'figure written' || echo 'no figure??')" \
        || { echo "  FAILED -- last 20 lines of the log:"; tail -20 "$T/log"; exit 1; }
fi
echo

echo "=== ready ==="
echo "  Enable image generation in a pipeline run with:"
echo
echo "    cd $REPO"
echo "    qsub -v VIZ=true env/gadi_run.sh"
echo "    # (combine with BAKTA=true,... as needed)"
echo
echo "  Or render one figure by hand:"
echo
echo "    module load singularity"
echo "    singularity exec -B /g/data/rg42,/scratch/rg42 $SIF \\"
echo "      python3 $REPO/scripts/pangraph_blocks.py \\"
echo "        --graph GRAPH_JSON --outdir OUTDIR"
