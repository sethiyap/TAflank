#!/bin/bash
#PBS -P rg42
#PBS -q copyq
#PBS -l walltime=00:30:00
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l jobfs=1GB
#PBS -l storage=gdata/rg42+scratch/rg42
#PBS -l wd
#PBS -N pangraph_setup
#PBS -j oe
#
# One-time setup for BUILD_PANGRAPH_GRAPH on Gadi: download the pangraph
# binary. Much simpler than env/gadi_bakta_setup.sh or env/gadi_viz_setup.sh
# -- pangraph 1.4.0 is a single statically-linked Rust binary (no Julia, no
# container, no database), so this is a plain file download, not a
# container pull or a multi-GB database fetch. Still has to happen
# somewhere with network, though -- copyq or a login node, never a compute
# queue.
#
#   cd /g/data/rg42/TAflank
#   qsub env/gadi_pangraph_setup.sh          # via copyq (recommended)
#   bash env/gadi_pangraph_setup.sh          # or directly on a login node
#
# Idempotent: re-running skips the download if the binary is already there
# and its version checks out.
#
# WHICH BUILD. The musl target (not gnu) is used deliberately: it is fully
# statically linked, so it runs on whatever glibc Gadi's compute nodes
# happen to have without needing it matched at build time -- the same
# "no runtime surprises" reasoning the native arm64 BLAST build on the
# laptop profile already follows.
set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
PGROOT="${PGROOT:-/g/data/rg42/pangraph}"
VERSION="${VERSION:-1.4.0}"
ARCH="${ARCH:-x86_64-unknown-linux-musl}"
URL="https://github.com/neherlab/pangraph/releases/download/${VERSION}/pangraph-${ARCH}"
BIN="$PGROOT/bin/pangraph"

mkdir -p "$PGROOT/bin"

echo "version : $VERSION"
echo "arch    : $ARCH"
echo "install : $BIN"
echo

echo "=== 1. download ==="
if [ -x "$BIN" ] && "$BIN" --version 2>/dev/null | grep -q "$VERSION"; then
    echo "  already present: $BIN ($("$BIN" --version))"
else
    curl -fsSL -o "$BIN" "$URL"
    chmod +x "$BIN"
    echo "  downloaded: $BIN ($(du -h "$BIN" | cut -f1))"
fi
echo

echo "=== 2. version check ==="
"$BIN" --version
echo

echo "=== 3. smoke test on one real locus set, if available ==="
LOCI_DIR=$(ls -d "$REPO"/loci/plasmids "$REPO"/run_all*/loci 2>/dev/null | head -1 || true)
if [ -z "$LOCI_DIR" ]; then
    echo "  skipped: no loci dir found yet (run the pipeline first)"
else
    T="${PBS_JOBFS:-/tmp}/pangraph_smoke"
    rm -rf "$T"; mkdir -p "$T"
    cat "$LOCI_DIR"/*.fa > "$T/combined.fa"
    "$BIN" build "$T/combined.fa" -o "$T/graph.json" > "$T/log" 2>&1 \
        && echo "  OK  $(/opt/anaconda3/bin/python3.12 -c "import json; print(len(json.load(open('$T/graph.json'))['blocks']), 'blocks')" 2>/dev/null || echo 'graph.json written')" \
        || { echo "  FAILED -- last 20 lines of the log:"; tail -20 "$T/log"; exit 1; }
fi
echo

echo "=== ready ==="
echo "  Enable pangraph in a pipeline run with:"
echo
echo "    cd $REPO"
echo "    qsub -v PANGRAPH=true env/gadi_run.sh"
echo "    # (combine with BAKTA=true,VIZ=true as needed)"
echo
echo "  Or build one graph by hand:"
echo
echo "    $BIN build LOCI_DIR/*.fa -o OUTPUT/graph.json"
