#!/bin/bash
#PBS -P rg42
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l storage=scratch/rg42
#PBS -l wd
#PBS -N fetch_blastdb
#PBS -j oe
#
# Download NCBI pre-formatted BLAST databases on Gadi.
#
# WHY A JOB AND NOT A LOGIN-NODE ONE-LINER
#   nt_prok is 128 GB compressed (42 volumes, 8.96 M sequences, 503 Gbp) and
#   nt_viruses another 68.7 GB. That is hours of transfer. Login nodes are
#   resource-limited and sessions drop. copyq is Gadi's data-transfer queue and
#   is one of the few with outbound network access -- standard compute queues
#   have none.
#
# WHY NOT $HOME
#   Gadi enforces a 10 GiB per-user home quota that `df` does not show. A
#   download there dies with curl error 23, "Failed writing body", which looks
#   like a network fault but is the write being refused. /scratch has petabytes.
#
# CAVEAT -- /scratch IS PURGED. Files not accessed for ~100 days are deleted.
#   Either re-touch the database periodically, or accept re-downloading it.
#   /g/data is persistent but neither ei36 (65 GB free) nor rg42 (54 GB free)
#   can hold nt_prok.
#
# Usage:  qsub env/gadi_fetch_blastdb.sh
#   override the project or target:
#         qsub -P ei36 -v DEST=/scratch/ei36/$USER/blastdb env/gadi_fetch_blastdb.sh

set -euo pipefail
PROJ="${PROJ:-rg42}"
DEST="${DEST:-/scratch/${PROJ}/${USER}/blastdb}"
# taxdb is small and gives staxids/sscinames in the BLAST output, which lets the
# pipeline classify records with NO network call -- essential on compute nodes.
DBS="${DBS:-nt_prok taxdb}"

mkdir -p "$DEST"
cd "$DEST"
echo "destination : $DEST"
echo "databases   : $DBS"
df -h "$DEST" | tail -1

module load blast+ 2>/dev/null || echo "note: no blast+ module; using update_blastdb.pl from PATH"

# --num_threads 8 parallelises the decompression. If you see repeated
# "(23) Failed writing body", check quota FIRST -- it is almost never the network.
update_blastdb.pl --decompress --num_threads 8 --source gcp $DBS

echo "=== done ==="
du -sh "$DEST"
ls "$DEST"/*.nal "$DEST"/*.nin 2>/dev/null | head
echo
echo "Add to your environment (and to the pipeline's --blastdb):"
echo "  export BLASTDB=$DEST"
echo "  ./tools/nf run main.nf -profile gadi --blastdb $DEST/nt_prok"
