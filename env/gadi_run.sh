#!/bin/bash
#PBS -P rg42
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l jobfs=20GB
#PBS -l storage=gdata/rg42+scratch/rg42
#PBS -l wd
#PBS -N taflank
#PBS -j oe
#
# Run the whole TAflank pipeline inside one PBS job.
#
# Everything runs with Nextflow's `local` executor inside this job, so there is
# no head process to keep alive on a login node. For the nt_prok sweep use
# -profile gadi instead, where tasks fan out across nodes.
#
# Everything is set with `qsub -v NAME=VALUE`. Several at once are comma
# separated, with NO spaces anywhere in the -v argument.
#
# Usage:
#   qsub env/gadi_run.sh                              # defaults below
#   qsub -v REPLICON=chromosome env/gadi_run.sh       # chromosomes instead
#   qsub -v TIER=all,REPLICON=all env/gadi_run.sh     # keep everything
#
# WHAT TO SEARCH (the subject; pick one, not both):
#   qsub -v GENOMES=/path/to/fasta_dir env/gadi_run.sh      # build a db from FASTAs
#   qsub -v DB=/scratch/rg42/blastdb/nt_prok env/gadi_run.sh # pre-formatted db
#
# WHAT TO SEARCH WITH (the query / bait):
#   qsub -v BAIT=/g/data/rg42/TAflank/anchors/relBE_seed.fa env/gadi_run.sh
#
# WHERE THE OUTPUT GOES (default /scratch/rg42/$USER/taflank_run):
#   qsub -v OUT=/scratch/rg42/$USER/run2 env/gadi_run.sh
#   Tables land in $OUT/results, windows in $OUT/loci, Nextflow work in
#   $PBS_JOBFS (discarded when the job ends). Keep OUT on /scratch, not in
#   $HOME -- home is capped at 10 GB. /scratch is purged after ~100 days of no
#   access, so copy anything worth keeping to /g/data/rg42/.
#
# OTHER SETTINGS:
#   FLANK=20000   window size each side of the bait
#   REPO=...      a different checkout of TAflank
#   BAKTA=true BAKTA_DB=/g/data/rg42/bakta_db/db-light
#                 full annotation. Run env/gadi_bakta_setup.sh first: it pulls
#                 the container AND downloads the database, both of which need
#                 network. Compute nodes have none, so neither can happen here.
#   BAKTA_SIF=... override the container image. Default: the newest
#                 bakta-*.sif beside BAKTA_DB, which is where the installer
#                 puts it. Nextflow is given this FILE, never the docker://
#                 URI -- handed a URI it tries to pull from quay.io on the
#                 compute node and dies with "network is unreachable".
#
#   A different bait needs ENDS changed with it:
#   qsub -v BAIT=/path/parDE.fa,ENDS=IS26 env/gadi_run.sh
#
#   ENDS names the IS elements that bracket the core, so it decides
#   COMPLETE / DECAYED / HALF-CORE / NO CORE. ISEsa1/ISEsa2 are specific to the
#   relBE metal-resistance lineage; left in place for another bait, every locus
#   reports NO CORE and nothing warns you.
#
#   SUBREGIONS is OPTIONAL and cannot change what is filtered: FULL and
#   PARALOGUE are decided before it is consulted, qcov>=95 already excludes
#   single-gene hits, and tier_a=FULL drops ORPHAN_* and PARTIAL alike. It only
#   adds cov_<part> reporting columns, and only earns its keep with TIER=all,
#   where a half hit resolves to ORPHAN_<part> instead of a flat PARTIAL.
#   SUBREGIONS=none is a perfectly good way to run this.

set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
OUT="${OUT:-/scratch/rg42/$USER/taflank_run}"
DB="${DB:-}"                      # empty => build from GENOMES
GENOMES="${GENOMES:-$REPO/data/public_16/fasta}"
BAIT="${BAIT:-$REPO/anchors/relBE_seed.fa}"
SUBREGIONS="${SUBREGIONS-relB:1-240,relE:240-527}"   # no :- , '' is valid
ENDS="${ENDS:-ISEsa1,ISEsa2}"
TIER="${TIER:-FULL}"
REPLICON="${REPLICON:-plasmid}"
FLANK="${FLANK:-20000}"
BAKTA="${BAKTA:-false}"
BAKTA_DB="${BAKTA_DB:-/g/data/rg42/bakta_db/db-light}"
BAKTA_SIF="${BAKTA_SIF:-}"

# Fail now, with a readable message, rather than 20 minutes into the queue.
[ -s "$BAIT" ] || { echo "BAIT not found or empty: $BAIT" >&2; exit 2; }
if [ -n "$DB" ]; then
    [ -s "$DB.nin" ] || [ -s "$DB.nal" ] || {
        echo "DB does not look like a BLAST nucleotide database: $DB" >&2
        echo "give the path WITHOUT an extension, e.g. /scratch/rg42/blastdb/nt_prok" >&2
        exit 2; }
else
    [ -d "$GENOMES" ] || { echo "GENOMES is not a directory: $GENOMES" >&2; exit 2; }
fi

if [ "$BAKTA" = "true" ]; then
    [ -s "$BAKTA_DB/version.json" ] || {
        echo "BAKTA_DB does not look like a Bakta database: $BAKTA_DB" >&2
        echo "run env/gadi_bakta_setup.sh first (needs network: copyq or login)" >&2
        exit 2; }
    if [ -z "$BAKTA_SIF" ]; then
        # newest bakta-*.sif sitting beside the database, where the installer
        # leaves it; `ls -t` so a re-pull at a newer tag wins automatically
        # `|| true`: with `set -e` a failed glob would exit here silently,
        # before the readable message below ever ran
        BAKTA_SIF=$(ls -t "$(dirname "$BAKTA_DB")"/bakta-*.sif 2>/dev/null | head -1) || true
    fi
    [ -s "$BAKTA_SIF" ] || {
        echo "no Bakta container image found next to $BAKTA_DB" >&2
        echo "compute nodes have NO outbound network, so Nextflow cannot pull" >&2
        echo "one. Run env/gadi_bakta_setup.sh, or pass BAKTA_SIF=/path.sif" >&2
        exit 2; }
fi

module purge
module load nextflow/24.04.5
module load blast/2.11.0
[ "$BAKTA" = "true" ] && module load singularity

export NXF_OPTS='-Xms512m -Xmx4g'
export NXF_HOME="${NXF_HOME:-/scratch/rg42/$USER/.nextflow}"
export TMPDIR="${PBS_JOBFS:-/scratch/rg42/$USER/tmp}"
mkdir -p "$OUT" "$TMPDIR" "$NXF_HOME"

cd "$REPO"
echo "repo     : $REPO   ($(git log --oneline -1))"
echo "outdir   : $OUT"
echo "database : ${DB:-<built from $GENOMES>}"
echo "bait     : $BAIT  ($(grep -c '^>' "$BAIT") sequence(s))"
echo "subregion: ${SUBREGIONS:-<none>}"
echo "ends     : $ENDS"
echo "filters  : tier_a=$TIER  replicon_type=$REPLICON  flank=$FLANK"
echo "bakta    : $BAKTA${BAKTA:+  ${BAKTA_SIF:+image=$BAKTA_SIF}}"
echo

ARGS=(-profile gadi_single
      --outdir   "$OUT/results"
      --loci_dir "$OUT/loci"
      --bait     "$BAIT"
      --subregions "$SUBREGIONS"
      --vehicle_ends "$ENDS"
      --tier_a   "$TIER"
      --replicon_type "$REPLICON"
      --flank    "$FLANK"
      -work-dir  "$TMPDIR/work"
      -resume)

if [ -n "$DB" ]; then
    # Pre-formatted database (nt_prok): no FASTAs ship with it, so sequence is
    # pulled back out of the database itself with blastdbcmd.
    ARGS+=(--blastdb "$DB" --genome_dir '' --blast_task megablast --taxonomy)
else
    ARGS+=(--genome_dir "$GENOMES" --gbk_dir "$REPO/data/public_16/gbk")
fi

[ "$BAKTA" = "true" ] && ARGS+=(--run_bakta true --bakta_db "$BAKTA_DB" \
                                --bakta_container "$BAKTA_SIF")

nextflow run main.nf "${ARGS[@]}"

echo
echo "=== results ==="
ls -la "$OUT/results"/*.tsv 2>/dev/null | awk '{print "  ", $5, $9}'
echo
[ -f "$OUT/results/summary.tsv" ] && column -t -s$'\t' "$OUT/results/summary.tsv" | cut -c1-120
