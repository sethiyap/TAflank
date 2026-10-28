#!/bin/bash
# One-time TAflank setup on Gadi. Run on a LOGIN NODE -- it needs outbound
# network, which standard compute queues do not have.
#
#   bash env/gadi_setup.sh
#
# Fetches what the repository deliberately does not carry: the genome records
# (252 MB, recoverable by accession) and the two reference databases. Idempotent
# -- re-running skips anything already present.
set -euo pipefail

REPO="${REPO:-/g/data/rg42/TAflank}"
cd "$REPO"
module load blast/2.11.0 2>/dev/null || true
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
KEY=""; [ -f ~/.ncbi_api_key ] && KEY="&api_key=$(cat ~/.ncbi_api_key)"
DELAY=$([ -n "$KEY" ] && echo 0.12 || echo 0.35)

echo "=== 1. ISfinder BLAST database ==="
if [ -s db/isfinder/ISfinder.nin ]; then
    echo "  already built"
else
    # NO -parse_seqids: the community mirror contains duplicate ids (e.g.
    # IS1004_IS200/IS605_IS200), which makeblastdb rejects outright. We only
    # ever BLAST against this database, never retrieve from it by id, so
    # sequence ids are not needed -- and this matches how the local copy was
    # built.
    makeblastdb -in db/isfinder/IS.fna -dbtype nucl \
        -out db/isfinder/ISfinder -title "ISfinder mirror" > /dev/null
    echo "  built from IS.fna ($(grep -c '^>' db/isfinder/IS.fna) sequences)"
fi

echo "=== 2. PlasmidFinder database ==="
if [ -s db/plasmidfinder/plasmidfinder.fa ]; then
    echo "  already present"
else
    mkdir -p db/plasmidfinder
    rm -rf db/plasmidfinder/_src
    git clone -q --depth 1 \
        https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git \
        db/plasmidfinder/_src
    python3 - <<'PY'
import glob, os
with open('db/plasmidfinder/plasmidfinder.fa', 'w') as out:
    n = 0
    for f in sorted(glob.glob('db/plasmidfinder/_src/*.fsa')):
        grp = os.path.basename(f)[:-4]
        for line in open(f):
            if line.startswith('>'):
                n += 1
                out.write('>%s|%s\n' % (grp, line[1:].strip()))
            else:
                out.write(line)
print('  built: %d replicon sequences' % n)
PY
fi

echo "=== 3. genome records (16 accessions) ==="
mkdir -p data/public_16/fasta data/public_16/gbk
need=0
while IFS=$'\t' read -r acc rest; do
    case "$acc" in \#*|'') continue;; esac
    base="${acc%%.*}"
    for fmt in fasta gb; do
        ext=$([ "$fmt" = fasta ] && echo fna || echo gbk)
        dir=$([ "$fmt" = fasta ] && echo fasta || echo gbk)
        out="data/public_16/$dir/${acc}.$ext"
        [ -s "$out" ] && continue
        need=1
        curl -sS --retry 3 --retry-delay 2 \
            "${EU}?db=nuccore&id=${acc}&rettype=${fmt}&retmode=text${KEY}" -o "$out"
        if [ ! -s "$out" ] || ! head -1 "$out" | grep -qE '^(>|LOCUS)'; then
            echo "  FAILED $acc ($fmt)"; rm -f "$out"
        fi
        sleep "$DELAY"
    done
done < data/public_16/accessions.tsv
echo "  fasta: $(ls data/public_16/fasta/*.fna 2>/dev/null | wc -l | tr -d ' ')/16   gbk: $(ls data/public_16/gbk/*.gbk 2>/dev/null | wc -l | tr -d ' ')/16"
[ "$need" -eq 0 ] && echo "  (all already present)"

echo
echo "=== ready ==="
echo "  qsub env/gadi_run.sh"
