#!/bin/bash
#PBS -P rg42
#PBS -q copyq
#PBS -l walltime=05:00:00
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l jobfs=20GB
#PBS -l storage=gdata/rg42+scratch/rg42
#PBS -l wd
#PBS -N bakta_setup
#PBS -j oe
#
# One-time Bakta setup on Gadi: pull the container, download the database.
#
#   cd /g/data/rg42/TAflank
#   qsub env/gadi_bakta_setup.sh          # via copyq (recommended)
#   bash env/gadi_bakta_setup.sh          # or directly on a login node
#
# It MUST run somewhere with outbound network. Standard compute queues (normal,
# express) have none, so this is copyq -- which exists for exactly this -- or a
# login node. Login nodes cap CPU time at 1800 s (`ulimit -t`), and unpacking
# the archive is xz decompression, so copyq is the safer of the two.
#
# Idempotent: re-running skips whatever is already in place.
#
# If PBS rejects the resource request above (copyq limits are not identical on
# every project), do not fight it -- run the same script on a login node:
#   bash env/gadi_bakta_setup.sh
#
# Options:
#   qsub -v TYPE=full env/gadi_bakta_setup.sh       # see the quota note below
#   qsub -v DBROOT=/scratch/rg42/$USER/bakta_db env/gadi_bakta_setup.sh
#
# WHICH DATABASE. `light` is 1.3 GB compressed, 3.9 GB unpacked. `full` is
# 30 GB compressed, 84 GB unpacked -- which does NOT fit the rg42 /g/data
# allocation (100 GiB total, ~54 GiB free when this was written), so a full
# database has to go on /scratch and will be purged after ~100 days of no
# access.
#
# Start with `light`. No conclusion in this project depends on the choice: CDS
# coordinates come from Pyrodigal and are database-independent, the vehicle is
# named by ISfinder, cargo by the curated gene set, and gene families (Stage 3)
# by MMseqs2 sequence clustering. Only product NAMES differ -- `light` keeps all
# non-coding features but only the PSCC/UniRef50 layer for CDS, dropping the
# exact-match layers that carry gene symbols and allele detail.
#
# The cost is real but secondary: 153 of 278 annotation-named CDSs across the
# seven plasmid loci are already `hypothetical protein` from GenBank (55%,
# ranging 22% in CP006642 to 76% in MF344582), so `full` would put names on some
# of them. Worth it only if characterising unnamed regions becomes the goal --
# e.g. CP028197's 13.4 kb gap between hok_gef and the tra operon. Switching is a
# re-download, not a redesign: identical pipeline, identical output columns.
# See PIPELINE_DESCRIPTION.md, "light vs full -- measured, not assumed".
set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
DBROOT="${DBROOT:-/g/data/rg42/bakta_db}"
TYPE="${TYPE:-light}"
TAG="${TAG:-1.12.2--pyhdfd78af_0}"
SIF="$DBROOT/bakta-${TAG}.sif"
DBDIR="$DBROOT/$([ "$TYPE" = light ] && echo db-light || echo db)"

case "$TYPE" in light|full) ;; *) echo "TYPE must be light or full" >&2; exit 2;; esac

module purge 2>/dev/null || true
module load singularity
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-/scratch/rg42/$USER/singularity_cache}"
mkdir -p "$DBROOT" "$SINGULARITY_CACHEDIR"

echo "database type : $TYPE"
echo "container     : quay.io/biocontainers/bakta:$TAG"
echo "install into  : $DBROOT"
echo "singularity   : $(singularity --version)"
echo

echo "=== 0. space check ==="
lquota 2>/dev/null | awk 'NR<=3 || /rg42/' || echo "  (lquota unavailable)"
need=$([ "$TYPE" = light ] && echo 6 || echo 120)
echo "  this needs roughly ${need} GiB free at $DBROOT (download + unpacked)"
echo

echo "=== 1. container ==="
if [ -s "$SIF" ]; then
    echo "  already pulled: $SIF"
else
    singularity pull "$SIF" "docker://quay.io/biocontainers/bakta:$TAG"
    echo "  pulled: $SIF  ($(du -h "$SIF" | cut -f1))"
fi

# -B is required: Singularity does not mount /g/data or /scratch by default,
# and without them the container cannot see where the database is going.
BAKTA=(singularity exec -B /g/data/rg42,/scratch/rg42 "$SIF")
echo "  bakta version: $("${BAKTA[@]}" bakta --version 2>&1 | tail -1 || true)"
echo

echo "=== 2. database ($TYPE) ==="
if [ -s "$DBDIR/version.json" ]; then
    echo "  already present: $DBDIR"
else
    # Creates $DBROOT/db-light (light) or $DBROOT/db (full) -- the name is
    # fixed by bakta/db.py, not something we choose. It also sets up the
    # AMRFinderPlus database, so no separate amrfinder_update is needed.
    "${BAKTA[@]}" bakta_db download --output "$DBROOT" --type "$TYPE"
fi
echo

echo "=== 3. verify ==="
[ -s "$DBDIR/version.json" ] || { echo "  FAILED: no $DBDIR/version.json" >&2; exit 1; }
echo "  db path    : $DBDIR"
echo "  db version : $(cat "$DBDIR/version.json")"
echo "  db size    : $(du -sh "$DBDIR" | cut -f1)"
if [ -d "$DBDIR/amrfinderplus-db" ]; then
    echo "  amrfinderplus-db present"
else
    echo "  NOTE no amrfinderplus-db subdirectory -- if bakta later complains," \
         "run: ${BAKTA[*]} amrfinder_update --force_update --database $DBDIR/amrfinderplus-db/"
fi
echo

echo "=== 4. smoke test on one real window ==="
W=$(ls "$REPO"/loci/plasmids/*.fa 2>/dev/null | head -1 || true)
if [ -z "$W" ]; then
    echo "  skipped: no window FASTA in $REPO/loci/plasmids (run the pipeline first)"
else
    T="${PBS_JOBFS:-/tmp}/bakta_smoke"
    rm -rf "$T"; mkdir -p "$T"
    "${BAKTA[@]}" bakta --db "$DBDIR" --output "$T" --prefix smoke \
        --threads 1 --keep-contig-headers --force \
        --skip-plot --skip-crispr "$W" > "$T/log" 2>&1 \
        && echo "  OK  $(grep -c $'\tCDS\t' "$T/smoke.gff3" || true) CDS called in $(basename "$W")" \
        || { echo "  FAILED -- last 20 lines of the bakta log:"; tail -20 "$T/log"; exit 1; }
fi
echo

echo "=== ready ==="
echo "  Add the database to a pipeline run with:"
echo
echo "    cd $REPO"
echo "    qsub -v BAKTA=true,BAKTA_DB=$DBDIR env/gadi_run.sh"
echo
echo "  Or annotate one window by hand:"
echo
echo "    module load singularity"
echo "    singularity exec -B /g/data/rg42,/scratch/rg42 $SIF \\"
echo "      bakta --db $DBDIR --output OUTDIR --prefix NAME \\"
echo "            --threads 4 --keep-contig-headers --force \\"
echo "            --skip-plot --skip-crispr WINDOW.fa"
