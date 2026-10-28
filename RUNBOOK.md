# TAflank — step by step

Run one stage at a time, look at what it produced, then move on. Every step is
a standalone script with its own CLI, so nothing here needs Nextflow and no
output is hidden in a work directory.

**Verified:** running all eleven steps by hand produces results identical to
`./tools/nf run main.nf` — 45 IS calls, 7 loci, 17 replicon hits, same context
classes. The pipeline is a wrapper around exactly these commands.

```bash
export BLAST_BIN=/Users/pooja/opt/ncbi-blast-2.17.0+/bin   # local
# module load blast/2.11.0                                  # Gadi
export OUT=/tmp/step        # anywhere you like
mkdir -p $OUT
```

---

## 1 · Build the search database

```bash
./scripts/build_blastdb.sh data/public_16/fasta $OUT/public16 "my run"
```

→ `$OUT/public16.*` and a `.dbmeta.tsv` sidecar recording source, sha256,
sequence count and BLAST version.

**Look at:** `cat $OUT/public16.dbmeta.tsv`
*Skip this step entirely if you are using a pre-formatted database
(`nt_prok`) — just point step 2 at it.*

---

## 2 · Bait search — Tier A, sequence

```bash
python3 scripts/bait_search.py \
    --bait anchors/relBE_seed.fa --db $OUT/public16 \
    --subregions relB:1-240,relE:240-527 \
    --out $OUT/bait_hits.tsv \
    --raw-out $OUT/blast/01_raw.tsv --kept-out $OUT/blast/02_tier_a.tsv
```

```
  [1] raw                 15 HSPs
  [2] after Tier A        14 HSPs
  FULL 14 · PARALOGUE 1
```

**Look at:** `blast/01_raw.tsv` is the verbatim blastn output with the exact
command as a comment. `02_tier_a.tsv` is the same columns after sequence
filtering, so `diff` shows precisely what identity removed.

**Try changing:** `--full-pident 90` to loosen, `--paralogue-max-pident 70`
to move the superfamily cut.

---

## 3 · Record metadata — Tier B input

```bash
python3 scripts/fetch_metadata.py \
    --accessions $OUT/bait_hits.tsv --col sacc --out $OUT/record_metadata.tsv
```

→ `chromosome 7 · plasmid 7`

**Look at:** the `replicon_type_evidence` column — `esummary_genome` is a
structured call, `title_*` is inferred. Cached per accession, so re-running is
free; add `--offline` to rebuild from cache, or `--from-blast` to parse
`stitle` and never touch the network.

---

## 4 · Filter — Tier B, record

```bash
python3 scripts/select_loci.py \
    --hits $OUT/bait_hits.tsv --meta $OUT/record_metadata.tsv \
    --out $OUT/loci_selected.tsv --flank 20000 \
    --tier-a FULL --replicon-type plasmid \
    --summary-out $OUT/filter_summary.tsv \
    --blast-out $OUT/blast/03_selected.tsv
```

```
  15 hits -> 7 selected
  deselected CP006642  tier_a=PARALOGUE
  deselected CP010512  replicon_type=chromosome   ...
```

**Look at:** `filter_summary.tsv` for the cascade, `loci_selected.tsv` for
every hit with `selected` and `deselect_reason`. Nothing is deleted.

**Try changing:** `--replicon-type chromosome`, or `--tier-a all
--replicon-type all` to keep everything and classify later.

---

## 5 · Cut the windows

```bash
python3 scripts/extract_loci.py \
    --selected $OUT/loci_selected.tsv --fasta-dir data/public_16/fasta \
    --out-dir $OUT/loci --out-table $OUT/loci.tsv --flank 20000
```

→ 7 FASTAs of 40,527 bp, bait at 20001–20527, relB on `+` in all of them.

**Look at:** the FASTA headers — they carry `orig=`, `flank_up/dn=`,
`truncated_up/dn=`, `wrapped=`, `whole_replicon=`.

**Other sources:** `--blastdb <path>` pulls sequence out of a BLAST database
(required for `nt_prok`, which ships no FASTAs); `--efetch` fetches from NCBI
for imported web-BLAST hits.

---

## 6 · IS elements — ISfinder, no length cutoff

```bash
for f in $OUT/loci/*.fa; do
  python3 scripts/is_search.py --locus "$f" --db db/isfinder/ISfinder \
      --out "$OUT/is/$(basename "$f" .fa).is.tsv"
done
```

```
  CP018961...: 63 HSPs -> 9 tiled (3 fragment)
```

**Look at:** `cov_class` (complete / partial / fragment), `is_fragment`,
`offset_from_bait`, and `nested_in`. The 64 bp ISEsa1 fragment in CP018961 at
offset −6739 is the one to find — it is the evidence the whole ancestral-
vehicle test rests on.

---

## 7 · Cargo genes

```bash
for f in $OUT/loci/*.fa; do
  python3 scripts/cargo_search.py --locus "$f" --db assets/cargo_genes.fa \
      --out "$OUT/cargo/$(basename "$f" .fa).cargo.tsv"
done
```

**Look at:** `offset_from_bait`. `czc` sits at −1184 in every locus that
carries it; `hok_gef` at +596/+600.

---

## 8 · Replicon typing — PlasmidFinder

```bash
for a in $(awk -F'\t' 'NR>1{print $2}' $OUT/loci.tsv | sort -u); do
  python3 scripts/plasmidfinder.py \
      --replicon data/public_16/fasta/${a}.1.fna \
      --db db/plasmidfinder/plasmidfinder.fa --accession "$a" \
      --out "$OUT/pf/${a}.pf.tsv"
done
```

> Runs on the **whole parent replicon**, never the window — Inc *rep* genes
> are rarely within 20 kb of the bait.

---

## 9 · Context — Tier C

```bash
python3 scripts/classify_context.py \
    --loci $OUT/loci.tsv --is-dir $OUT/is --cargo-dir $OUT/cargo \
    --out $OUT/context.tsv
```

```
  COMPLETE 3 · DECAYED 1 · HALF-CORE 2 · NO CORE 1
  with the czc+cytb metal module: 4/7
```

**Look at:** `vehicle_interrupted` — CP018961 is COMPLETE by coverage but
split by ISKpn21, which is a different thing from intact.

---

## 10 · Feature track — the figure's input

```bash
for f in $OUT/loci/*.fa; do
  b=$(basename "$f" .fa); acc=$(echo "$b" | cut -d_ -f1)
  python3 scripts/gene_blocks.py --locus "$f" \
      --is "$OUT/is/$b.is.tsv" --cargo "$OUT/cargo/$b.cargo.tsv" \
      --gbk "data/public_16/gbk/${acc}.1.gbk" \
      --out "$OUT/blocks/$b.blocks.tsv"
done
```

**Look at:** `colour_key` vs `colour_source`. A CDS annotated
"ISNCY family transposase" is coloured **ISEsa1**, because annotation names IS
elements by family only. Add `--gff` once Bakta output exists.

---

## 11 · Collate and report

```bash
python3 scripts/collate.py --is-dir $OUT/is --pf-dir $OUT/pf \
    --loci $OUT/loci.tsv --out-mge $OUT/mge_hits.tsv \
    --out-rep $OUT/replicons.tsv --out-summary $OUT/summary.tsv

python3 scripts/run_report.py --results $OUT --out $OUT/run_report.txt
```

**Look at:** `run_report.txt` — the whole run on one page.

---

## When you are happy with the flow

```bash
./tools/nf run main.nf -profile local            # laptop
qsub env/gadi_run.sh                             # Gadi
```

Same commands, same results, just orchestrated — with `-resume` so an
unchanged step is never recomputed.
