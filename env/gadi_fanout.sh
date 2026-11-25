#!/bin/bash
#PBS -P rg42
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=4
#PBS -l mem=8GB
#PBS -l storage=gdata/rg42+scratch/rg42
#PBS -l wd
#PBS -N taflank_fanout
#PBS -j oe
#
# Submit the TAflank Nextflow HEAD process to normal, which submits and polls
# one PBS job per task via -profile gadi. The tasks themselves still land on
# normal/other queues per configs/gadi.config; only the orchestrator runs
# here.
#
# Was copyq (outbound network, but capped at ncpus=1). Switched to normal
# because this run needs no runtime network -- meta_source defaults to
# 'blast' (no NCBI calls) and the Bakta/viz/pangraph containers/binaries are
# all pre-staged by the one-time gadi_*_setup.sh scripts -- and normal lets
# the head job request enough CPUs (ncpus=4 above) for the SELECT_LOCI /
# COLLATE / FETCH_METADATA steps, which run under Nextflow's `local` executor
# INSIDE this job and inherit `cpus = 4` from configs/gadi.config. On copyq's
# ncpus=1 cap those steps crash immediately with "Process requirement
# exceeds available CPUs -- req: 4; avail: 1". If a future run sets
# meta_source=ncbi (live E-utilities calls) or needs to pull a container it
# hasn't staged yet, switch this back to copyq (and ncpus back to 1, or check
# whether copyq allows more).
#
# Was walltime=10:00:00, then 24:00:00: at ~1031 loci x 5 fanned-out
# process types, plus this account's existing queue backlog, the head job
# (which only submits and polls -- no heavy compute of its own, ~1GB RSS,
# a few CPU-minutes) kept hitting its own walltime with no error, not even
# at 24h -- that run got to IS_SEARCH/BAKTA/PLASMIDFINDER 100% cached and
# CARGO_SEARCH/GENE_BLOCKS 2-3 tasks short before being killed
# ("job killed: walltime 86408 exceeded limit 86400"). -resume means a
# re-submission after a walltime kill just continues (cached work is not
# re-run), so 48:00:00 (normal queue's max) gives headroom rather than
# needing yet another resubmit for the last few stragglers + viz steps.
#
# This is the nt_prok-scale sweep. For the small curated set use
# env/gadi_run.sh (-profile gadi_single, one PBS job, no fan-out) instead.
#
# Usage:
#   qsub -v REPLICON=plasmid    env/gadi_fanout.sh
#   qsub -v REPLICON=chromosome env/gadi_fanout.sh
#   qsub -v REPLICON=chromosome,ORGANISM="Escherichia coli,Klebsiella pneumoniae" env/gadi_fanout.sh
#
# Run env/gadi_bait_search.sh FIRST (also on copyq or a login node) to see
# how many hits the search space actually has before committing to the full
# fan-out -- see that script's own header for why.
#
# WHAT TO SEARCH (a pre-formatted BLAST db; this script does not build one
# from FASTAs -- use env/gadi_run.sh for that):
#   qsub -v DB=/scratch/rg42/blastdb/nt_prok env/gadi_fanout.sh   # default
#
# WHERE THE OUTPUT GOES (default /scratch/rg42/$USER/taflank_fanout_<REPLICON>):
#   qsub -v OUT=/scratch/rg42/$USER/run2 env/gadi_fanout.sh
#   Nextflow's -work-dir lives under $OUT/work -- on /scratch, NOT jobfs.
#   Fan-out tasks run as separate PBS jobs on other nodes and need a shared
#   filesystem to find each other's work directories; jobfs is local to
#   whichever single node this head job happens to land on and the other
#   tasks could never see it.
#
# ONE-TIME SETUP (each needs network -- copyq or a login node -- because
# compute nodes have none to pull a container or download a binary):
#   qsub env/gadi_blast_setup.sh
#   qsub env/gadi_bakta_setup.sh
#   qsub env/gadi_viz_setup.sh
#   qsub env/gadi_pangraph_setup.sh
# This script assumes all four have already been run once; the full-stack
# flags below (BAKTA/VIZ/PANGRAPH) default to true, and BLAST_BIN_DIR is
# always required, each failing fast with a readable message if its
# container/binary isn't found.
#
# WHY gadi_blast_setup.sh: blast/2.11.0 -- the newest blastn Gadi has as a
# MODULE -- cannot read nt_prok's current def-line format and dies with
# CSerialException::eOverflow. There is no newer module to switch to
# (blast/2.10.1 is the only other option, older still), so the fix is a
# container, via a wrapper bin directory (see that script's own header).
#
# OTHER SETTINGS (same meaning as env/gadi_run.sh -- see its header for the
# long version of each):
#   BAIT=...      FLANK=20000      REPO=...      TIER=FULL
#   ORGANISM=...  comma-separated, case-insensitive substrings matched
#                 against the organism field, e.g. "Escherichia coli,
#                 Klebsiella pneumoniae". Default 'all' keeps every organism --
#                 use this AFTER a replicon-type sweep finishes, to narrow the
#                 already-selected loci down to specific species, rather than
#                 re-running BAIT_SEARCH.
#   SUBREGIONS=... (default relB:1-240,relE:240-527)
#   ENDS=ISEsa1,ISEsa2
#   BLAST_BIN_DIR=/g/data/rg42/blast_env/bin (default)
#   BAKTA_DB=/g/data/rg42/bakta_db/db-light   BAKTA_SIF=...
#   VIZ_SIF=...   PANGRAPH_BIN=...
#   Set BAKTA=false / VIZ=false / PANGRAPH=false to drop any one of the three
#   full-stack pieces instead of the base tables.

set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
DB="${DB:-/scratch/rg42/blastdb/nt_prok}"
REPLICON="${REPLICON:-plasmid}"
ORGANISM="${ORGANISM:-all}"
OUT="${OUT:-/scratch/rg42/$USER/taflank_fanout_$REPLICON}"
BAIT="${BAIT:-$REPO/anchors/relBE_seed.fa}"
SUBREGIONS="${SUBREGIONS-relB:1-240,relE:240-527}"   # no :- , '' is valid
ENDS="${ENDS:-ISEsa1,ISEsa2}"
TIER="${TIER:-FULL}"
FLANK="${FLANK:-20000}"
BLAST_BIN_DIR="${BLAST_BIN_DIR:-/g/data/rg42/blast_env/bin}"
BAKTA="${BAKTA:-true}"
BAKTA_DB="${BAKTA_DB:-/g/data/rg42/bakta_db/db-light}"
BAKTA_SIF="${BAKTA_SIF:-}"
VIZ="${VIZ:-true}"
VIZ_SIF="${VIZ_SIF:-}"
PANGRAPH="${PANGRAPH:-true}"
PANGRAPH_BIN="${PANGRAPH_BIN:-}"

# Fail now, with a readable message, rather than hours into the fan-out.
[ -s "$BAIT" ] || { echo "BAIT not found or empty: $BAIT" >&2; exit 2; }
[ -s "$DB.nin" ] || [ -s "$DB.nal" ] || {
    echo "DB does not look like a BLAST nucleotide database: $DB" >&2
    echo "give the path WITHOUT an extension, e.g. /scratch/rg42/blastdb/nt_prok" >&2
    exit 2; }
[ -x "$BLAST_BIN_DIR/blastn" ] || {
    echo "no working blastn found in $BLAST_BIN_DIR" >&2
    echo "blast/2.11.0 (the Gadi module) cannot read $DB -- CSerialException::eOverflow" >&2
    echo "run env/gadi_blast_setup.sh first, or pass BLAST_BIN_DIR=/path" >&2
    exit 2; }

if [ "$BAKTA" = "true" ]; then
    [ -s "$BAKTA_DB/version.json" ] || {
        echo "BAKTA_DB does not look like a Bakta database: $BAKTA_DB" >&2
        echo "run env/gadi_bakta_setup.sh first (needs network: copyq or login)" >&2
        exit 2; }
    if [ -z "$BAKTA_SIF" ]; then
        BAKTA_SIF=$(ls -t "$(dirname "$BAKTA_DB")"/bakta-*.sif 2>/dev/null | head -1) || true
    fi
    [ -s "$BAKTA_SIF" ] || {
        echo "no Bakta container image found next to $BAKTA_DB" >&2
        echo "run env/gadi_bakta_setup.sh, or pass BAKTA_SIF=/path.sif" >&2
        exit 2; }
fi

if [ "$VIZ" = "true" ]; then
    if [ -z "$VIZ_SIF" ]; then
        VIZ_SIF=$(ls -t /g/data/rg42/viz_env/*.sif 2>/dev/null | head -1) || true
    fi
    [ -s "$VIZ_SIF" ] || {
        echo "no viz container image found in /g/data/rg42/viz_env" >&2
        echo "run env/gadi_viz_setup.sh, or pass VIZ_SIF=/path.sif" >&2
        exit 2; }
fi

if [ "$PANGRAPH" = "true" ]; then
    if [ -z "$PANGRAPH_BIN" ]; then
        PANGRAPH_BIN=$(dirname "$(ls -t /g/data/rg42/pangraph/bin/pangraph 2>/dev/null | head -1)") || true
    fi
    [ -x "$PANGRAPH_BIN/pangraph" ] || {
        echo "no pangraph binary found in /g/data/rg42/pangraph/bin" >&2
        echo "run env/gadi_pangraph_setup.sh, or pass PANGRAPH_BIN=/dir" >&2
        exit 2; }
fi

module purge
module load nextflow/24.04.5
module load singularity
# NOT `module load blast/2.11.0` -- that module can't read nt_prok (see
# BLAST_BIN_DIR check above). --blast_bin below overrides it completely.

export NXF_OPTS='-Xms512m -Xmx4g'
export NXF_HOME="${NXF_HOME:-/scratch/rg42/$USER/.nextflow}"
mkdir -p "$OUT" "$NXF_HOME"

cd "$REPO"
echo "repo     : $REPO   ($(git log --oneline -1))"
echo "outdir   : $OUT"
echo "database : $DB"
echo "bait     : $BAIT  ($(grep -c '^>' "$BAIT") sequence(s))"
echo "subregion: ${SUBREGIONS:-<none>}"
echo "ends     : $ENDS"
echo "filters  : tier_a=$TIER  replicon_type=$REPLICON  organism=$ORGANISM  flank=$FLANK"
echo "blast_bin: $BLAST_BIN_DIR"
echo "bakta    : $BAKTA  ${BAKTA_SIF:+image=$BAKTA_SIF}"
echo "viz      : $VIZ  ${VIZ_SIF:+image=$VIZ_SIF}"
echo "pangraph : $PANGRAPH  ${PANGRAPH_BIN:+bin=$PANGRAPH_BIN}"
echo

ARGS=(-profile gadi
      --blastdb  "$DB" --genome_dir '' --blast_task megablast --taxonomy
      --blast_bin "$BLAST_BIN_DIR"
      --outdir   "$OUT/results"
      --figdir   "$OUT/figures"
      --loci_dir "$OUT/loci"
      --bait     "$BAIT"
      --subregions "$SUBREGIONS"
      --vehicle_ends "$ENDS"
      --tier_a   "$TIER"
      --replicon_type "$REPLICON"
      --organism "$ORGANISM"
      --flank    "$FLANK"
      -work-dir  "$OUT/work"
      -resume)

[ "$BAKTA" = "true" ] && ARGS+=(--run_bakta true --bakta_db "$BAKTA_DB" \
                                --bakta_container "$BAKTA_SIF")
[ "$VIZ" = "true" ] && ARGS+=(--run_viz true --viz_container "$VIZ_SIF")
[ "$PANGRAPH" = "true" ] && ARGS+=(--run_pangraph true --pangraph_bin "$PANGRAPH_BIN")

nextflow run main.nf "${ARGS[@]}"

echo
echo "=== results ==="
ls -la "$OUT/results"/*.tsv 2>/dev/null | awk '{print "  ", $5, $9}'
echo
[ -f "$OUT/results/summary.tsv" ] && column -t -s$'\t' "$OUT/results/summary.tsv" | cut -c1-120
