#!/bin/bash
#PBS -P rg42
#PBS -q copyq
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/rg42+scratch/rg42
#PBS -l wd
#PBS -N blast_setup
#PBS -j oe
#
# One-time setup for BAIT_SEARCH against a pre-formatted external database
# (nt_prok): pull a modern blast+ container, because the newest blastn
# available as a Gadi MODULE (blast/2.11.0; blast/2.10.1 is the only other
# option, older still) cannot read nt_prok's current ASN.1 def-line format
# and dies with:
#   CSerialException::eOverflow ... error while reading seqid
#   ( at Blast-def-line-set.[].[].seqid.[].[].gi)
# -- a documented old-client/new-database mismatch, not a config problem.
# configs/local.config already pins blast_bin to a 2.17.0 install for this
# exact reason (`/Users/pooja/opt/ncbi-blast-2.17.0+/bin`); this script gets
# Gadi an equivalent version, just via container instead of a local install.
#
#   cd /g/data/rg42/TAflank
#   qsub env/gadi_blast_setup.sh             # via copyq (recommended)
#   bash env/gadi_blast_setup.sh             # or directly on a login node
#
# Idempotent: re-running skips whatever is already pulled.
#
# WHAT THIS BUILDS: not just a container, but a directory of wrapper SCRIPTS
# -- /g/data/rg42/blast_env/bin/{blastn,blastdbcmd,makeblastdb} -- each a
# one-line `exec singularity exec ... blast-X.Y.Z.sif <tool> "$@"` shim. This
# matters because every script in this repo that calls blastn/blastdbcmd
# resolves the binary via BLAST_BIN as a plain PATH-style directory prefix
# (scripts/toolpath.py's resolve(): `os.path.join(bin_dir, name)`), not as
# an arbitrary command -- it has no way to express "run this through
# singularity" itself. Pointing BLAST_BIN at this directory instead gets
# every blast call (BAIT_SEARCH's blastn, but also EXTRACT_LOCI's and
# PLASMIDFINDER's blastdbcmd against the same nt_prok -- same bug, same db)
# the modern binary, with zero changes to main.nf or any config.
set -euo pipefail

BLASTROOT="${BLASTROOT:-/g/data/rg42/blast_env}"
BLAST_IMAGE="${BLAST_IMAGE:-docker://ncbi/blast:2.17.0}"
DB="${DB:-/scratch/rg42/blastdb/nt_prok}"
BAIT="${BAIT:-/g/data/rg42/TAflank/anchors/relBE_seed.fa}"
TAG=$(basename "$BLAST_IMAGE" | tr ':' '-')
SIF="$BLASTROOT/${TAG}.sif"
BIN="$BLASTROOT/bin"

module purge 2>/dev/null || true
module load singularity
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-/scratch/rg42/$USER/singularity_cache}"
mkdir -p "$BLASTROOT" "$BIN" "$SINGULARITY_CACHEDIR"

echo "image         : $BLAST_IMAGE"
echo "install into  : $SIF"
echo "wrappers in   : $BIN"
echo "singularity   : $(singularity --version)"
echo

echo "=== 1. container ==="
if [ -s "$SIF" ]; then
    echo "  already pulled: $SIF"
else
    singularity pull "$SIF" "$BLAST_IMAGE"
    echo "  pulled: $SIF  ($(du -h "$SIF" | cut -f1))"
fi
echo

echo "=== 2. wrapper scripts ==="
for tool in blastn blastdbcmd makeblastdb; do
    cat > "$BIN/$tool" <<WRAP
#!/bin/bash
exec singularity exec -B /g/data/rg42,/scratch/rg42 "$SIF" $tool "\$@"
WRAP
    chmod +x "$BIN/$tool"
    echo "  wrote $BIN/$tool"
done
echo

echo "=== 3. version check ==="
"$BIN/blastn" -version
echo

echo "=== 4. smoke test: the exact command that failed, through the new binary ==="
if [ -s "$DB.nin" ] || [ -s "$DB.nal" ]; then
    T="${PBS_JOBFS:-/tmp}/blast_smoke"
    rm -rf "$T"; mkdir -p "$T"
    BLAST_BIN="$BIN" python3 /g/data/rg42/TAflank/scripts/bait_search.py \
        --bait "$BAIT" --db "$DB" --out "$T/bait_hits.tsv" \
        --subregions 'relB:1-240,relE:240-527' \
        --full-qcov 95 --paralogue-max-pident 90 \
        --task megablast --evalue 1e-5 --max-target-seqs 100000 \
        --num-threads 1 --taxonomy \
        --raw-out "$T/01_raw.tsv" --kept-out "$T/02_tier_a.tsv" \
        > "$T/log" 2>&1 \
        && echo "  OK  $(wc -l < "$T/bait_hits.tsv") hit row(s) -- no CSerialException" \
        || { echo "  FAILED -- last 30 lines of the log:"; tail -30 "$T/log"; exit 1; }
else
    echo "  skipped: $DB does not look like a BLAST db (set DB=... to point at nt_prok)"
fi
echo

echo "=== ready ==="
echo "  Point BLAST_BIN at the wrapper directory in a pipeline run with:"
echo
echo "    qsub -v BLAST_BIN_DIR=$BIN,... env/gadi_fanout.sh   # already the default"
echo
echo "  Or by hand:"
echo
echo "    BLAST_BIN=$BIN python3 scripts/bait_search.py ..."
