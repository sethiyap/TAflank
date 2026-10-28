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
# Usage:
#   qsub env/gadi_run.sh                              # defaults below
#   qsub -v REPLICON=chromosome env/gadi_run.sh       # chromosomes instead
#   qsub -v TIER=all,REPLICON=all env/gadi_run.sh     # keep everything
#   qsub -v DB=/scratch/rg42/blastdb/nt_prok env/gadi_run.sh   # sweep database

set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
OUT="${OUT:-/scratch/rg42/$USER/taflank_run}"
DB="${DB:-}"                      # empty => build from GENOMES
GENOMES="${GENOMES:-$REPO/data/public_16/fasta}"
TIER="${TIER:-FULL}"
REPLICON="${REPLICON:-plasmid}"
FLANK="${FLANK:-20000}"
BAKTA="${BAKTA:-false}"
BAKTA_DB="${BAKTA_DB:-/g/data/rg42/bakta_db/db-light}"

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
echo "filters  : tier_a=$TIER  replicon_type=$REPLICON  flank=$FLANK"
echo "bakta    : $BAKTA"
echo

ARGS=(-profile gadi_single
      --outdir   "$OUT/results"
      --loci_dir "$OUT/loci"
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
