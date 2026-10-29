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

### The three things you can change

Everything below is written out in full, so you can see exactly where each
piece comes from and where it goes. There are only three decisions:

| | what it is | the default used below |
|---|---|---|
| **the bait** | the sequence you search *with* — the seed | `/g/data/rg42/TAflank/anchors/relBE_seed.fa` (527 bp: relB + relE) |
| **the database** | the sequences you search *in* | built in STEP 1 from the 16 records, into `/scratch/rg42/ps1744/myrun/mydatabase` |
| **the output folder** | where everything lands | `/scratch/rg42/ps1744/myrun` |

To use your own, substitute the path wherever it appears. Nothing else changes:

```bash
# your own bait -- any FASTA, one or more sequences
--bait /scratch/rg42/ps1744/my_own_seed.fa

# a database you built yourself, or a downloaded one like nt_prok
# (give the path with NO file extension)
--db /scratch/rg42/ps1744/blastdb/nt_prok

# a different output folder -- make it first
mkdir -p /scratch/rg42/ps1744/run2
```

Put the output on **`/scratch`**, not in your home directory: home is capped at
10 GB and a full run will blow through it. Note that `/scratch` is purged after
about 100 days of no access, so copy anything you want to keep to
`/g/data/rg42/`.

For the one-command version of all of this, see
[Telling it what to search, and what to search with](#telling-it-what-to-search-and-what-to-search-with)
at the end — the same three choices become `qsub -v BAIT=...,DB=...,OUT=...`.

### What the output files mean

`OUTPUTS.md` in the repo describes every table and every column:

```bash
less /g/data/rg42/TAflank/OUTPUTS.md
```

### Looking at a table

Every "look at it" block below uses the same small helper, which picks columns
**by name**:

```bash
python3 /g/data/rg42/TAflank/scripts/tsvcols.py <file> <column> <column> ...
```

Run it with no column names and it lists what the file contains:

```bash
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/cargo/CP028197__CP028197__1__plus.cargo.tsv
```

Ask for a name that is not there and it says so, and lists the real ones. It
replaces `cut -f2,8,9,10 | column -t`, which looks equivalent but silently
hands back *different* columns as soon as a new one is added to a table — which
is exactly what happened when `clipped_at_window_edge` arrived at position 10
of the cargo table and `offset_from_bait` quietly dropped out of the view. A
name cannot go wrong that way.

A glob matching several files gives one combined table with the header printed
once, so this works too:

```bash
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/cargo/*.cargo.tsv \
  locus_id family offset_from_bait
```

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
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/record_metadata.tsv \
  accession organism replicon_type replicon_name replicon_type_evidence
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
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/is/CP018961*.is.tsv \
  element family win_start win_end query_cov cov_class is_fragment \
  clip_cause clipped_missing_bp offset_from_bait
```

**The one to find:** a row saying `ISEsa1 ... 64 ... fragment ... -6739`.
That 64-base scrap is the single most important piece of evidence in the whole
project — it is what shows the same vehicle carried relBE into each plasmid.

**Check the `clip_cause` column before believing any `fragment`.** A short
alignment has two possible explanations and they mean opposite things:

| `clip_cause` | what happened | what it means |
|---|---|---|
| `-` | the element really is a short remnant in open DNA | **evidence** — this is what Stage 5 counts |
| `window` | the element is whole, but our ±20 kb cut sliced it | **not evidence** — we just stopped looking there |
| `replicon` | the replicon itself ends there | real truncation, a signal in its own right |

The ISEsa1 row at −6739 says `-`, so it is genuine. Two rows in the whole set
say `window` (in CP028197 and CP028537) and both are unrelated IS110-family
elements sitting on the window boundary. The script prints the split for you:

```
... (1 fragment, of which 1 cut by the window edge; 0 usable as remnants)
```

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
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/cargo/CP028197*.cargo.tsv \
  family pident cov clip_cause clipped_missing_bp offset_from_bait
```

**What to notice:** `czc` sits at position **−1184** in every plasmid that has
it. Not approximately — exactly. That is a strong hint they share an ancestor.

The reading order upstream → downstream is `cytb` (−1960), `czc` (−1184),
`relB` (−2), `relE` (+241), `hok_gef` (+596) — the complete core, with a second
toxin–antitoxin system right behind it — then `tra` conjugation genes at +14 kb.
`tra` genes are **not** transposable elements: they move the whole plasmid into
another cell, they do not cut and paste themselves into DNA. They are backbone
context. The vehicle that inserted relBE is in the IS table from step 6.

Note the last `tra` row: `cov` 84.1% but `clip_cause` = `window`, 499 bp
missing. The gene is not degraded — it runs off the end of the window.

**Now the column that carries the actual result**, which the view above hides:

```bash
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/cargo/CP028197*.cargo.tsv \
  family source cov offset_from_bait
```

`source` names the record each reference gene came from. For CP028197 the
`cytb`, `czc`, `relB` and `relE` matches are 100% identical to copies taken
from **four different records across three genera** — *Citrobacter*,
*Enterobacter* and *Escherichia*, two chromosomes and two plasmids. Four
adjacent genes do not stay 100% identical across three genera by descent. The
whole ~2.5 kb block moved recently, horizontally, as one unit.

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
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/context.tsv \
  locus_id vehicle_class vehicle_ends vehicle_interrupted \
  vehicle_edge_clipped n_is_fragment n_is_fragment_real cargo_families
```

| class | meaning |
|---|---|
| COMPLETE | both ends of the vehicle are there |
| DECAYED | both ends there but one is broken |
| HALF-CORE | only one end left |
| NO CORE | the vehicle is gone — relBE arrived some other way |

---

# STEP 9b — Annotate the windows with Bakta  *(optional, and not yet runnable)*

**This is where Bakta sits: between classifying the loci and building the
drawing input.** It is optional, and right now it cannot run on Gadi. Both
reasons are below, with what to do about them.

Bakta and ISfinder are **complementary, not alternatives**:

| | what it gives you | what it cannot do |
|---|---|---|
| **ISfinder** (STEP 6) | the IS **element name** — `ISEsa1`, by homology to a named reference | nothing outside its database; no genes |
| **Bakta** (this step) | every CDS in the window, uniformly called | name the vehicle. It labels IS elements by **family only** — "ISNCY family transposase" — which cannot distinguish ISEsa1 from ISKpn21 |

So the vehicle call comes from ISfinder and the gene track from Bakta, and
STEP 10 merges them with the element name winning on ≥50% overlap.

### What you get without it

STEP 10 already works: it falls back to each record's **own GenBank
annotation**, and `annot_source` in the output reads `genbank`. That is why the
blocks tables exist today. It is good enough for 15 of the 16 records.

The exception is **LT991956**, which has **0 CDS features** in its GenBank
record. Check for yourself:

```bash
python3 /g/data/rg42/TAflank/scripts/tsvcols.py \
  /scratch/rg42/ps1744/myrun/blocks/LT991956*.blocks.tsv \
  feature_type name annot_source
```

Every row says `IS` or `cargo`; there are no `CDS` rows. That locus cannot
contribute a gene track to the figure until it is annotated.

### Why it will not run yet

Two things are missing on Gadi, and the pipeline is already wired for both:

1. **The database is not there.** `/g/data/rg42/bakta_db/` does not exist. The
   light database is about 1.5 GB and must be downloaded from a node with
   network access — a login node or `copyq`, never a compute node.
2. **The container has to be pulled once.** It is now named in the config
   (`quay.io/biocontainers/bakta:1.12.2--pyhdfd78af_0`, checked against the
   real quay.io tag list), and Gadi has `module load singularity`. But the pull
   also needs network, so it cannot happen inside a compute job.

Bakta cannot run on your laptop at all, for an unrelated reason:
`amrfinderplus`, which Bakta hard-requires, has **no macOS build** in bioconda
or conda-forge — not arm64, not x86_64, not noarch. Rosetta does not help.
That is why the GenBank fallback exists.

### To set it up (one time, from a login node)

```bash
mkdir -p /g/data/rg42/bakta_db
cd /g/data/rg42/bakta_db
module load singularity
export SINGULARITY_CACHEDIR=/scratch/rg42/$USER/singularity_cache
singularity pull --name bakta.sif \
  docker://quay.io/biocontainers/bakta:1.12.2--pyhdfd78af_0
singularity exec bakta.sif bakta_db download --output . --type light
```

Then run the pipeline with it switched on:

```bash
cd /g/data/rg42/TAflank
qsub -v BAKTA=true,BAKTA_DB=/g/data/rg42/bakta_db/db-light env/gadi_run.sh
```

I have not run any of this — the download and the pull are yours to start, and
the `db-light` folder name should be confirmed against what `bakta_db download`
actually creates before the `-v` path above is trusted.

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

This is the table the coloured-block figure will be drawn from. The
`annot_source` column records where the CDS calls came from: `genbank` with the
command above, `bakta` if you ran STEP 9b and swap `--gbk` for
`--gff /scratch/rg42/ps1744/myrun/bakta/$b.gff3`.

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

## Telling it what to search, and what to search with

Nothing is edited. Everything is set with `qsub -v NAME=VALUE`, and several at
once are comma separated with **no spaces anywhere** in the `-v` argument.

**What to search** (the subject). Pick one of these two, not both:

```bash
# (a) build a database from a directory of FASTA files
qsub -v GENOMES=/g/data/rg42/TAflank/data/public_16/fasta env/gadi_run.sh

# (b) use a database that is already formatted — give the path with NO extension
qsub -v DB=/scratch/rg42/ps1744/blastdb/nt_prok env/gadi_run.sh
```

`GENOMES` is the default, pointing at the 16 records. `DB` overrides it: with a
pre-formatted database there are no FASTA files to read, so each window is
pulled back out of the database itself with `blastdbcmd`.

**What to search with** (the query, i.e. the bait):

```bash
qsub -v BAIT=/g/data/rg42/TAflank/anchors/relBE_seed.fa env/gadi_run.sh
```

**Where the output goes.** The default is
`/scratch/rg42/ps1744/taflank_run` — note that this is *not* the
`myrun` folder the eleven steps above use, so a `qsub` run will not overwrite
what you built by hand.

```bash
qsub -v OUT=/scratch/rg42/ps1744/run2 env/gadi_run.sh
```

| goes where | what |
|---|---|
| `$OUT/results/` | every table in `OUTPUTS.md` |
| `$OUT/results/blast/` | the BLAST output before and after filtering |
| `$OUT/results/blocks/`, `is/`, `cargo/`, `pf/` | the per-locus tables |
| `$OUT/loci/` | the ±20 kb window FASTAs |
| `$PBS_JOBFS` | Nextflow's scratch — **discarded when the job ends** |

Keep `OUT` on `/scratch`, never in `$HOME` (capped at 10 GB). `/scratch` is
purged after about 100 days without access, so copy results worth keeping to
`/g/data/rg42/`.

**Everything else:** `FLANK=20000` sets the window size, `REPO=` points at a
different checkout, and `BAKTA=true` turns on annotation (see STEP 9b — it
needs setting up first). The full list:

```bash
sed -n '19,76p' /g/data/rg42/TAflank/env/gadi_run.sh
```

Both together:

```bash
qsub -v BAIT=/g/data/rg42/TAflank/anchors/relBE_seed.fa,DB=/scratch/rg42/ps1744/blastdb/nt_prok \
     env/gadi_run.sh
```

### If you change the bait, change `ENDS` with it

`ENDS` names the IS elements that bracket the core, and that is what
COMPLETE / DECAYED / HALF-CORE / NO CORE is measured against. `ISEsa1,ISEsa2`
belong to the **relBE metal-resistance lineage**. Leave them in place for a
different bait and every locus comes back NO CORE — not because the vehicle is
absent, but because you asked about the wrong vehicle, and nothing warns you.

It is load-bearing, not cosmetic. On the current seven loci the default gives
3 COMPLETE / 1 DECAYED / 2 HALF-CORE / 1 NO CORE; `ENDS=IS26,IS903B` on the
very same data gives 1 COMPLETE / 3 HALF-CORE / 3 NO CORE.

```bash
qsub -v BAIT=/g/data/rg42/TAflank/anchors/parDE.fa,ENDS=IS26 env/gadi_run.sh
```

### `SUBREGIONS` is optional — you do not have to supply coordinates

It cannot change what is filtered, for three independent reasons:

1. `FULL` and `PARALOGUE` are decided **before** subregions are consulted, so
   subregions can never alter either call;
2. `qcov >= 95` already excludes any single-gene hit anyway — relB alone is
   240 bp of a 527 bp bait, i.e. 45% coverage;
3. `tier_a=FULL` then drops `ORPHAN_*` and `PARTIAL` **identically**, so naming
   them differently changes nothing downstream.

What it adds is two reporting columns, `cov_relB` and `cov_relE`. On the one
non-FULL hit in the 16 records they read `relB 0.8%, relE 90.6%` — a bare
relE-family toxin with no antitoxin partner, which is the RelE/ParE superfamily
background confirmed rather than assumed. Useful to read; never decisive.

Omit it, or set `SUBREGIONS=none`, and nothing breaks:

```bash
qsub -v BAIT=/g/data/rg42/TAflank/anchors/parDE.fa,SUBREGIONS=none,ENDS=IS26 \
     env/gadi_run.sh
```

Name the parts only if you intend to run `TIER=all` and study half hits, where
the label becomes `ORPHAN_PARD` / `ORPHAN_PARE` instead of a flat `PARTIAL`:

```bash
qsub -v BAIT=/path/parDE.fa,SUBREGIONS=parD:1-250,parE:250-550,TIER=all,ENDS=IS26 \
     env/gadi_run.sh
```

### See all the settings

```bash
sed -n '19,76p' /g/data/rg42/TAflank/env/gadi_run.sh
```

The job prints every one of them at the top of its log before doing any work,
and stops immediately with a readable message if `BAIT` does not exist or `DB`
is not a BLAST database — rather than failing twenty minutes into the queue.
