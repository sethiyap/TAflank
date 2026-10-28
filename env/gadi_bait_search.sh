#!/bin/bash
#PBS -P rg42
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/rg42+scratch/rg42
#PBS -l wd
#PBS -N taflank_bait
#PBS -j oe
#
# TAflank -- Tier A scoping run: how many relBE copies are actually in nt_prok?
#
# WHY A JOB
#   Gadi login nodes enforce `ulimit -t 1800` -- 30 minutes of CPU per process.
#   A bait search against nt_prok (8,962,499 sequences, 502.7 Gbp) exceeds that
#   comfortably and gets killed, which looks like a dropped connection.
#
# WHY SCOPE FIRST
#   Every FULL hit becomes a window extraction, an IS search and a PlasmidFinder
#   run. On 16 genomes that was 7 loci and 20 tasks. If nt_prok yields 20,000
#   full-length hits that is 60,000 tasks. Find out N before committing.
#
# WHY megablast AND NOT -task blastn
#   Only FULL hits (>=95% identity) are carried forward, and megablast finds
#   those. Word size 11 exists to catch IS copies at 85-92% identity, which is
#   the IS search's problem, not the bait search's. Honest cost: the PARALOGUE
#   class is not reported at all rather than being logged -- a deviation from
#   "rejections are logged, never silent" that belongs in the methods.
#
# Usage:  qsub env/gadi_bait_search.sh
#   override:  qsub -v BAIT=anchors/other.fa,TASK=blastn env/gadi_bait_search.sh

set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
BLASTDB_DIR="${BLASTDB_DIR:-/scratch/rg42/blastdb}"
DB="${DB:-nt_prok}"
BAIT="${BAIT:-$REPO/anchors/relBE_seed.fa}"
TASK="${TASK:-megablast}"
OUT="${OUT:-/scratch/rg42/$USER/taflank_scope}"
SUBREG="${SUBREG:-relB:1-240,relE:240-527}"

module purge
module load blast/2.11.0

export BLASTDB="$BLASTDB_DIR"
mkdir -p "$OUT"
cd "$REPO"

echo "repo     : $REPO"
echo "database : $BLASTDB_DIR/$DB"
echo "bait     : $BAIT  ($(grep -v '^>' "$BAIT" | tr -d '\n' | wc -c | tr -d ' ') bp)"
echo "task     : $TASK   threads: ${PBS_NCPUS:-8}"
blastdbcmd -db "$DB" -info 2>/dev/null | head -3
echo

/usr/bin/time -v python3 scripts/bait_search.py \
    --bait "$BAIT" \
    --db "$DB" \
    --out "$OUT/bait_hits.tsv" \
    --subregions "$SUBREG" \
    --task "$TASK" \
    --num-threads "${PBS_NCPUS:-8}" \
    --max-target-seqs 5000000 \
    --taxonomy 2>&1 | grep -E "HSPs|FULL|PARALOGUE|ORPHAN|PARTIAL|Elapsed|Maximum resident" || true

echo
echo "=== scoping summary ==="
H="$OUT/bait_hits.tsv"
echo "  total HSPs          : $(( $(wc -l < "$H") - 1 ))"
echo "  distinct accessions : $(awk -F'\t' 'NR>1{print $2}' "$H" | sort -u | wc -l)"
echo
echo "  by Tier A class:"
awk -F'\t' 'NR>1{c[$(NF-1)]++} END{for(k in c) printf "    %-18s %d\n", k, c[k]}' "$H"
echo
echo "  FULL hits by replicon word in the subject title:"
awk -F'\t' 'NR>1 && $(NF-1)=="FULL"{
    t=tolower($15)
    if (t ~ /plasmid/) p++
    else if (t ~ /whole genome shotgun|scaffold|contig/) d++
    else if (t ~ /chromosome|complete genome|complete sequence/) c++
    else u++
} END{printf "    plasmid    %d\n    chromosome %d\n    draft      %d\n    unknown    %d\n", p,c,d,u}' "$H"
echo
echo "  written to: $H"
echo "  NEXT: review the numbers, set --tier-a / --replicon-type / --min-flank,"
echo "        then run the full pipeline with -resume so this search is not repeated."
