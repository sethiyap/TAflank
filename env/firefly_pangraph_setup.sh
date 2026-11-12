#!/usr/bin/env bash
# One-time setup for BUILD_PANGRAPH_GRAPH on Firefly: download the pangraph
# binary. Plain script, NOT a sbatch/qsub job -- unlike env/gadi_run.sh's
# copyq, Firefly has no confirmed network-capable queue yet (that is
# exactly what env/firefly_probe.sh's "outbound internet" section checks).
# Run this directly on whatever node firefly_probe.sh showed has outbound
# access (http=200 to github.com, specifically -- not just quay.io/NCBI).
#
#   bash env/firefly_pangraph_setup.sh
#
# If no node has outbound network, this cannot run anywhere on Firefly:
# download the binary elsewhere and scp it to PGROOT/bin/pangraph instead,
# then skip straight to the version-check step below.
#
# Idempotent: re-running skips the download if the binary is already there
# and its version checks out. Same musl-static reasoning as
# env/gadi_pangraph_setup.sh -- see that script's header.
set -euo pipefail

PGROOT="${PGROOT:-$HOME/opt/pangraph}"
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
    curl -fsSL -o "$BIN" "$URL" || {
        echo "  FAILED to download from $URL" >&2
        echo "  if this node has no outbound network, download it" >&2
        echo "  elsewhere and scp to $BIN, then rerun this script" >&2
        exit 1; }
    chmod +x "$BIN"
    echo "  downloaded: $BIN ($(du -h "$BIN" | cut -f1))"
fi
echo

echo "=== 2. version check ==="
"$BIN" --version
echo

echo "=== ready ==="
echo "  This binary needs to be visible from WHEREVER Nextflow tasks"
echo "  actually run on Firefly -- confirm \$PGROOT is on a shared"
echo "  filesystem (not node-local \$HOME), same caveat as this cluster's"
echo "  apptainer autoMounts in general (configs/firefly.config)."
echo
echo "  Enable pangraph in a pipeline run with:"
echo
echo "    nextflow run main.nf -profile firefly --run_pangraph true \\"
echo "      --pangraph_bin $PGROOT/bin ..."
