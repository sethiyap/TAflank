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
#   A different bait almost always needs SUBREGIONS and ENDS changed too:
#   qsub -v BAIT=/path/parDE.fa,SUBREGIONS=parD:1-250,ENDS=IS26 env/gadi_run.sh
#
#   SUBREGIONS names the parts of the bait, so per-gene coverage is reported
#   (use SUBREGIONS=none for a bait with no named parts). ENDS names the IS
#   elements that bracket the core; ISEsa1/ISEsa2 belong to the relBE
#   metal-resistance lineage, so leaving them in place for another bait makes
#   every locus report NO CORE.

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
echo "bakta    : $BAKTA"
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

[ "$BAKTA" = "true" ] && ARGS+=(--run_bakta true --bakta_db "$BAKTA_DB")

nextflow run main.nf "${ARGS[@]}"

echo
echo "=== results ==="
ls -la "$OUT/results"/*.tsv 2>/dev/null | awk '{print "  ", $5, $9}'
echo
[ -f "$OUT/results/summary.tsv" ] && column -t -s$'\t' "$OUT/results/summary.tsv" | cut -c1-120
