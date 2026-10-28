# TAflank — running it step by step on Gadi

Copy and paste one block at a time. Look at what came out. Move on.

Every path below is written out in full — there is nothing to set up, no
variables to define, and no editing needed.

---

## Before you start (once per login)

```bash
ssh gadi
module load blast/2.11.0
mkdir -p /scratch/rg42/ps1744/myrun
```

That is all. The databases and genome files are already in place.

---

# STEP 1 — Build the search database

```bash
/g/data/rg42/TAflank/scripts/build_blastdb.sh \
  /g/data/rg42/TAflank/data/public_16/fasta \
  /scratch/rg42/ps1744/myrun/mydatabase \
  "my first run"
```

**You should see:** `built /scratch/rg42/ps1744/myrun/mydatabase`

**Look at what it made:**

```bash
cat /scratch/rg42/ps1744/myrun/mydatabase.dbmeta.tsv
```

16 sequences, 46,893,952 bases. This is the haystack we will search.

---

# STEP 2 — Search for relBE

```bash
python3 /g/data/rg42/TAflank/scripts/bait_search.py \
  --bait /g/data/rg42/TAflank/anchors/relBE_seed.fa \
  --db /scratch/rg42/ps1744/myrun/mydatabase \
  --subregions relB:1-240,relE:240-527 \
  --out /scratch/rg42/ps1744/myrun/bait_hits.tsv \
  --raw-out /scratch/rg42/ps1744/myrun/01_raw.tsv \
  --kept-out /scratch/rg42/ps1744/myrun/02_after_sequence_filter.tsv
```

**You should see:**

```
  [1] raw                 15 HSPs
  [2] after Tier A        14 HSPs
15 HSPs
  FULL               14
  PARALOGUE          1
```

**What happened:** 15 matches found. 14 are real relBE. 1 is a look-alike from
the same protein family (67% identical) — not the gene we want.

**Look at the raw BLAST output:**

```bash
column -t -s$'\t' /scratch/rg42/ps1744/myrun/01_raw.tsv | cut -c1-110
```

**See exactly which match was removed:**

```bash
diff /scratch/rg42/ps1744/myrun/01_raw.tsv \
     /scratch/rg42/ps1744/myrun/02_after_sequence_filter.tsv
```

---

# STEP 3 — Find out what each match sits on

```bash
python3 /g/data/rg42/TAflank/scripts/fetch_metadata.py \
  --accessions /scratch/rg42/ps1744/myrun/bait_hits.tsv \
  --col sacc \
  --out /scratch/rg42/ps1744/myrun/record_metadata.tsv
```

**You should see:** `chromosome 7` and `plasmid 7`

**Look at it:**

```bash
cut -f1,17,18,19 /scratch/rg42/ps1744/myrun/record_metadata.tsv | column -t -s$'\t'
```

This asked NCBI what each record is. The `replicon_type_evidence` column says
how it knew — `esummary_genome` means NCBI stated it outright, `title_*` means
we read it off the description.

---

# STEP 4 — Keep only plasmids

```bash
python3 /g/data/rg42/TAflank/scripts/select_loci.py \
  --hits /scratch/rg42/ps1744/myrun/bait_hits.tsv \
  --meta /scratch/rg42/ps1744/myrun/record_metadata.tsv \
  --out /scratch/rg42/ps1744/myrun/loci_selected.tsv \
  --flank 20000 \
  --tier-a FULL \
  --replicon-type plasmid \
  --summary-out /scratch/rg42/ps1744/myrun/filter_summary.tsv \
  --blast-out /scratch/rg42/ps1744/myrun/03_after_plasmid_filter.tsv
```

**You should see:** `15 hits -> 7 selected`, then the reason each one was
dropped.

**Look at the whole filtering story on one screen:**

```bash
column -t -s$'\t' /scratch/rg42/ps1744/myrun/filter_summary.tsv
```

> Nothing is ever deleted. All 15 are still in `loci_selected.tsv`, each with a
> `selected` column saying yes or no and a `deselect_reason` saying why.

**To keep chromosomes instead**, change `--replicon-type plasmid` to
`--replicon-type chromosome`. **To keep everything**, use `--tier-a all
--replicon-type all`.

---

# STEP 5 — Cut out 20 kb either side of each hit

```bash
python3 /g/data/rg42/TAflank/scripts/extract_loci.py \
  --selected /scratch/rg42/ps1744/myrun/loci_selected.tsv \
  --fasta-dir /g/data/rg42/TAflank/data/public_16/fasta \
  --out-dir /scratch/rg42/ps1744/myrun/loci \
  --out-table /scratch/rg42/ps1744/myrun/loci.tsv \
  --flank 20000
```

**You should see:** `7 windows`

**Look at the files it made:**

```bash
ls /scratch/rg42/ps1744/myrun/loci/
head -1 /scratch/rg42/ps1744/myrun/loci/*.fa
```

Each is 40,527 bases long, with relBE sitting in the middle at position
20001–20527, and always pointing the same way.

---

# STEP 6 — Find the jumping genes (IS elements)

```bash
for f in /scratch/rg42/ps1744/myrun/loci/*.fa; do
  python3 /g/data/rg42/TAflank/scripts/is_search.py \
    --locus "$f" \
    --db /g/data/rg42/TAflank/db/isfinder/ISfinder \
    --out "/scratch/rg42/ps1744/myrun/is/$(basename "$f" .fa).is.tsv"
done
```

*(This runs the same command once for each of the 7 files.)*

**Look at one locus:**

```bash
cut -f2,3,5,6,10,11,12,13 \
  /scratch/rg42/ps1744/myrun/is/CP018961*.is.tsv | column -t -s$'\t'
```

**The one to find:** a row saying `ISEsa1 ... 64 ... fragment ... -6739`.
That 64-base scrap is the single most important piece of evidence in the whole
project — it is what shows the same vehicle carried relBE into each plasmid.

---

# STEP 7 — Find the cargo genes

```bash
for f in /scratch/rg42/ps1744/myrun/loci/*.fa; do
  python3 /g/data/rg42/TAflank/scripts/cargo_search.py \
    --locus "$f" \
    --db /g/data/rg42/TAflank/assets/cargo_genes.fa \
    --out "/scratch/rg42/ps1744/myrun/cargo/$(basename "$f" .fa).cargo.tsv"
done
```

**Look at one:**

```bash
cut -f2,8,9,10 \
  /scratch/rg42/ps1744/myrun/cargo/CP028197*.cargo.tsv | column -t -s$'\t'
```

**What to notice:** `czc` sits at position **−1184** in every plasmid that has
it. Not approximately — exactly. That is a strong hint they share an ancestor.

---

# STEP 8 — Work out the plasmid type

```bash
for a in CP006642 CP018961 CP023893 CP028197 CP028537 LT991956 MF344582; do
  python3 /g/data/rg42/TAflank/scripts/plasmidfinder.py \
    --replicon /g/data/rg42/TAflank/data/public_16/fasta/${a}.1.fna \
    --db /g/data/rg42/TAflank/db/plasmidfinder/plasmidfinder.fa \
    --accession ${a} \
    --out /scratch/rg42/ps1744/myrun/pf/${a}.pf.tsv
done
```

**You should see** one line per plasmid, e.g.
`CP028197: 3 replicon(s): IncHI2A,IncHI2,IncC`

---

# STEP 9 — Classify each locus

```bash
python3 /g/data/rg42/TAflank/scripts/classify_context.py \
  --loci /scratch/rg42/ps1744/myrun/loci.tsv \
  --is-dir /scratch/rg42/ps1744/myrun/is \
  --cargo-dir /scratch/rg42/ps1744/myrun/cargo \
  --out /scratch/rg42/ps1744/myrun/context.tsv
```

**You should see:** `COMPLETE 3 · DECAYED 1 · HALF-CORE 2 · NO CORE 1`

**Look at it:**

```bash
cut -f1,4,5,6,9 /scratch/rg42/ps1744/myrun/context.tsv | column -t -s$'\t'
```

| class | meaning |
|---|---|
| COMPLETE | both ends of the vehicle are there |
| DECAYED | both ends there but one is broken |
| HALF-CORE | only one end left |
| NO CORE | the vehicle is gone — relBE arrived some other way |

---

# STEP 10 — Build the drawing input

```bash
for f in /scratch/rg42/ps1744/myrun/loci/*.fa; do
  b=$(basename "$f" .fa)
  acc=$(echo "$b" | cut -d_ -f1)
  python3 /g/data/rg42/TAflank/scripts/gene_blocks.py \
    --locus "$f" \
    --is "/scratch/rg42/ps1744/myrun/is/$b.is.tsv" \
    --cargo "/scratch/rg42/ps1744/myrun/cargo/$b.cargo.tsv" \
    --gbk "/g/data/rg42/TAflank/data/public_16/gbk/${acc}.1.gbk" \
    --out "/scratch/rg42/ps1744/myrun/blocks/$b.blocks.tsv"
done
```

This is the table the coloured-block figure will be drawn from.

---

# STEP 11 — Put it all together

```bash
python3 /g/data/rg42/TAflank/scripts/collate.py \
  --is-dir /scratch/rg42/ps1744/myrun/is \
  --pf-dir /scratch/rg42/ps1744/myrun/pf \
  --loci /scratch/rg42/ps1744/myrun/loci.tsv \
  --out-mge /scratch/rg42/ps1744/myrun/mge_hits.tsv \
  --out-rep /scratch/rg42/ps1744/myrun/replicons.tsv \
  --out-summary /scratch/rg42/ps1744/myrun/summary.tsv
```

```bash
python3 /g/data/rg42/TAflank/scripts/run_report.py \
  --results /scratch/rg42/ps1744/myrun \
  --out /scratch/rg42/ps1744/myrun/run_report.txt
```

**The final answer, on one page:**

```bash
cat /scratch/rg42/ps1744/myrun/run_report.txt
```

---

# Once you are happy with all of this

The same eleven steps, run for you in one go:

```bash
cd /g/data/rg42/TAflank
qsub env/gadi_run.sh
```

It produces identical results — verified. The only difference is that you do
not watch each step go past.
