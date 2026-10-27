#!/usr/bin/env bash
# Build a nucleotide BLAST db from a directory of FASTAs, with a provenance
# sidecar. Replaces the /tmp-built databases whose sources no longer exist.
#
#   build_blastdb.sh <fasta_dir> <out_prefix> <db_title>
set -euo pipefail
SRC_DIR="$1"; OUT="$2"; TITLE="$3"
BLAST_BIN="${BLAST_BIN:-/Users/pooja/opt/ncbi-blast-2.17.0+/bin}"
mkdir -p "$(dirname "$OUT")" "$(dirname "$OUT")/_src"
CAT="${OUT}.source.fa"

: > "$CAT"
n=0
for f in "$SRC_DIR"/*.fna "$SRC_DIR"/*.fa "$SRC_DIR"/*.fasta; do
  [ -e "$f" ] || continue
  cat "$f" >> "$CAT"; n=$((n+1))
done
[ "$n" -gt 0 ] || { echo "no FASTA found in $SRC_DIR" >&2; exit 1; }

"$BLAST_BIN/makeblastdb" -in "$CAT" -dbtype nucl -out "$OUT" -title "$TITLE" \
  -parse_seqids -hash_index > /dev/null

{
  echo -e "key\tvalue"
  echo -e "db_title\t$TITLE"
  echo -e "source_dir\t$SRC_DIR"
  echo -e "source_fasta\t$CAT"
  echo -e "n_input_files\t$n"
  echo -e "n_sequences\t$(grep -c '^>' "$CAT")"
  echo -e "total_bp\t$(grep -v '^>' "$CAT" | tr -d '\n' | wc -c | tr -d ' ')"
  echo -e "sha256\t$(shasum -a 256 "$CAT" | cut -d' ' -f1)"
  echo -e "blast_version\t$("$BLAST_BIN/makeblastdb" -version | head -1 | awk '{print $2}')"
  echo -e "built_utc\t$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo -e "built_by\t$(whoami)@$(hostname -s)"
} > "${OUT}.dbmeta.tsv"

echo "built $OUT"
column -t -s$'\t' "${OUT}.dbmeta.tsv"
