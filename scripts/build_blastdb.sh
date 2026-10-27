#!/usr/bin/env bash
# Build a nucleotide BLAST db from a directory of FASTAs, with a provenance
# sidecar. Replaces the /tmp-built databases whose sources no longer exist.
#
#   build_blastdb.sh <fasta_dir> <out_prefix> <db_title>
#
# BLAST_BIN may name the directory holding makeblastdb. If it is empty the
# binary is resolved from PATH, which is what happens on a cluster where it
# comes from a container or a module.
set -euo pipefail
SRC_DIR="$1"; OUT="$2"; TITLE="$3"
BLAST_BIN="${BLAST_BIN:-}"

if [ -n "$BLAST_BIN" ] && [ -x "$BLAST_BIN/makeblastdb" ]; then
  MAKEDB="$BLAST_BIN/makeblastdb"
else
  MAKEDB="$(command -v makeblastdb || true)"
fi
[ -n "$MAKEDB" ] && [ -x "$MAKEDB" ] || {
  echo "build_blastdb.sh: makeblastdb not found." >&2
  echo "  set BLAST_BIN=/path/to/ncbi-blast-X/bin, or put it on PATH" >&2
  exit 1; }

mkdir -p "$(dirname "$OUT")"
CAT="${OUT}.source.fa"
: > "$CAT"
n=0
for f in "$SRC_DIR"/*.fna "$SRC_DIR"/*.fa "$SRC_DIR"/*.fasta; do
  [ -e "$f" ] || continue
  cat "$f" >> "$CAT"; n=$((n+1))
done
[ "$n" -gt 0 ] || { echo "no FASTA found in $SRC_DIR" >&2; exit 1; }

"$MAKEDB" -in "$CAT" -dbtype nucl -out "$OUT" -title "$TITLE" \
  -parse_seqids -hash_index > /dev/null

if command -v shasum >/dev/null 2>&1; then SHA="$(shasum -a 256 "$CAT" | cut -d' ' -f1)"
elif command -v sha256sum >/dev/null 2>&1; then SHA="$(sha256sum "$CAT" | cut -d' ' -f1)"
else SHA="unavailable"; fi

{
  echo -e "key\tvalue"
  echo -e "db_title\t$TITLE"
  echo -e "source_dir\t$SRC_DIR"
  echo -e "source_fasta\t$CAT"
  echo -e "n_input_files\t$n"
  echo -e "n_sequences\t$(grep -c '^>' "$CAT")"
  echo -e "total_bp\t$(grep -v '^>' "$CAT" | tr -d '\n' | wc -c | tr -d ' ')"
  echo -e "sha256\t$SHA"
  echo -e "makeblastdb\t$MAKEDB"
  echo -e "blast_version\t$("$MAKEDB" -version | head -1 | awk '{print $2}')"
  echo -e "built_utc\t$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo -e "built_by\t$(whoami)@$(hostname -s)"
} > "${OUT}.dbmeta.tsv"

echo "built $OUT"
