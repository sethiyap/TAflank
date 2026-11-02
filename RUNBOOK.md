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

# STEP 9b — Annotate the windows with Bakta  *(optional)*

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

Two things are missing on Gadi. The pipeline is already wired for both, so
installing them is all that is left:

1. **The database is not there** — `/g/data/rg42/bakta_db/` does not exist.
2. **The container has not been pulled.** It is named in the config as
   `quay.io/biocontainers/bakta:1.12.2--pyhdfd78af_0`, read off the real
   quay.io tag list, and Gadi has `module load singularity` (3.11.3).

Both steps need **outbound network**, which standard compute queues do not
have. So they go through `copyq`, or a login node.

### Installing it — one command

```bash
cd /g/data/rg42/TAflank
qsub env/gadi_bakta_setup.sh
```

That script pulls the container, downloads the database, verifies it, and
annotates one of your real windows as a smoke test before declaring success.
It is idempotent, so re-running it skips whatever is already in place. If PBS
rejects the resource request, run the same script on a login node instead:

```bash
bash env/gadi_bakta_setup.sh
```

Watch it with `qstat -u $USER`. When it finishes it prints the exact `qsub`
line to use next.

### Why the installer has to pull the container, and the pipeline must not

Gadi's compute queues have **no outbound network**. The installer runs on
`copyq`, which does, and leaves two things on disk:

```
/g/data/rg42/bakta_db/db-light                       the database
/g/data/rg42/bakta_db/bakta-1.12.2--pyhdfd78af_0.sif the container image
```

If Nextflow is handed `docker://quay.io/biocontainers/bakta:...` it will try
to fetch it itself, on a compute node, and the run dies part-way through:

```
Failed to pull singularity image
  pinging container registry quay.io: dial tcp ...:443: connect: network is unreachable
```

That failure is not confined to the BAKTA step. Nextflow kills the tasks
already in flight, so IS_SEARCH, CARGO_SEARCH and PLASMIDFINDER are left
half-finished and the output directory looks plausible but is incomplete —
`results/` with no `context.tsv` and no `blocks/` is the signature.

So `env/gadi_run.sh` resolves the **local `.sif` file** and passes that as
`--bakta_container`. It looks for the newest `bakta-*.sif` sitting beside
`BAKTA_DB`, which is exactly where the installer puts it, and refuses to
submit if it cannot find one. Override with `BAKTA_SIF=/path/to/image.sif`
if you keep it somewhere else.

### Which database — start with `light`

| | compressed | unpacked | fits in `/g/data/rg42`? |
|---|---|---|---|
| **`light`** | 1.3 GB | **3.9 GB** | yes — about 54 GiB free |
| `full` | 30 GB | **84 GB** | **no** — the whole rg42 allocation is 100 GiB |

**Nothing the project concludes depends on which one you pick:**

| what you need | comes from | changes with `full`? |
|---|---|---|
| CDS **coordinates** | Pyrodigal — gene calling, independent of the database | no |
| the **vehicle** (ISEsa1/ISEsa2) | ISfinder homology, STEP 6 | no |
| **cargo** (czc, cytb, relB, relE, hok_gef, tra) | the curated 333-sequence set, STEP 7 | no |
| **gene families** for the figure | MMseqs2 sequence clustering (Stage 3, not built) | no |
| gene **product names** | the Bakta database | **yes** |

Only the last row moves. And uniform *coordinates* were the point of deciding
to re-annotate at all — GenBank gives anywhere from 0 to 5,765 features
depending on who deposited the record. `light` fixes that completely.

### What `full` would actually buy

`light` keeps every non-coding feature but, for CDSs, only the **PSCC** layer
(UniRef50 clusters). It drops UPS/IPS/PSC — the exact-match layers that supply
gene symbols, allele-level detail and `WP_*`/UniRef100 cross-references — so
product names come out coarser.

That is worth weighing, because the current baseline is already weak. Measured
on the blocks tables built from GenBank annotation:

| locus | CDS named by annotation | "hypothetical protein" | |
|---|---|---|---|
| CP006642 | 32 | 7 | 22% |
| CP018961 | 39 | 19 | 49% |
| CP023893 | 60 | 27 | 45% |
| CP028197 | 51 | 33 | 65% |
| CP028537 | 50 | 32 | 64% |
| LT991956 | 0 | 0 | — (no CDS at all) |
| MF344582 | 46 | 35 | 76% |
| **total** | **278** | **153** | **55%** |

The 153 counts every key containing "hypothetical": 151 are exactly
`hypothetical protein` and 2 are `orf; hypothetical protein`. Reproduce it with:

```bash
awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
  $h["feature_type"]=="CDS" && $h["colour_source"]=="annotation" {print $h["colour_key"]}' \
  /scratch/rg42/ps1744/myrun/blocks/*.blocks.tsv > /tmp/keys.txt

wc -l < /tmp/keys.txt              # 278 annotation-named CDS
grep -ci hypothetical /tmp/keys.txt  # 153 of them unnamed
sort /tmp/keys.txt | uniq -c | sort -rn | head
```

That awk reads the header and looks columns up **by name**, so it keeps working
when columns are added — and unlike `tsvcols.py` piped into `cut`, it does not
break on product names that contain spaces.

So **55% of annotated CDSs are already unnamed**, and `full` would put names on
some of them. That matters if you want to know what sits in CP028197's unnamed
13.4 kb gap between `hok_gef` and the `tra` operon. It does not matter for the
vehicle, the cargo, or the offsets.

Switching later is a **re-download, not a redesign** — same script, same
pipeline, same output columns, only better product strings:

```bash
qsub -v TYPE=full,DBROOT=/scratch/rg42/$USER/bakta_db env/gadi_bakta_setup.sh
```

(`full` has to go on `/scratch`, which has room but is purged after ~100 days.)

### What neither database fixes

`gene_blocks.py` uses the product string as the colour key, so all 153
"hypothetical protein" CDSs collapse onto **one colour**. A better database
shrinks that pile; it cannot empty it.

The real fix is the Stage 3 MMseqs2 clustering — colour by **sequence** family,
not by whatever name the annotator happened to emit. An 84 GB download is not a
substitute for it; they solve different problems.

The download also sets up the AMRFinderPlus database itself, so there is no
separate `amrfinder_update` step.

### Then run the pipeline with it on

```bash
cd /g/data/rg42/TAflank
qsub -v BAKTA=true,BAKTA_DB=/g/data/rg42/bakta_db/db-light env/gadi_run.sh
```

`db-light` is not a name I chose — `bakta_db download --output DIR --type light`
always creates `DIR/db-light` (and `DIR/db` for `full`); the name is fixed in
Bakta's own `db.py`.

### Why it cannot be installed on your laptop

For a different reason entirely: `amrfinderplus`, which Bakta hard-requires,
has **no macOS build at all** — not in bioconda, not in conda-forge, not
arm64, not x86_64, not noarch. Rosetta does not help, because there is no
x86_64 build either. The only local route would be Docker, and that daemon is
off. This is why the GenBank fallback exists and why the merge in STEP 10 could
be built and tested before Bakta was available anywhere.

**I have not run the installation** — the container pull and the database
download are yours to start.

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

# STEP 11b — Cluster the proteins into gene families

Without this the CDS track is one flat grey: `gene_blocks.py` keys colour on
the product string, and 127 of the Bakta-named CDS are `hypothetical protein`.

```bash
python3 scripts/gene_families.py \
  --gff   run_bakta/bakta \
  --loci  run_bakta/loci \
  --out   results/gene_families.tsv \
  --check
```

`--check` asserts that relB and relE translate to the sequences recorded in
CLAUDE.md (79 aa `MGSINLRIDDEL…TLDEL`, 95 aa `MAYFLDFDERAL…VKRIL`). Keep it on
— it is a free proof that the proteins are right before anything is built on
them.

**Proteins are translated from the GFF3, not read from Bakta's `.faa`.** The
BAKTA process originally declared only the GFF3 as output, so publishDir never
copied the proteins and they died with the work directory (`$PBS_JOBFS`, gone
at job end). `main.nf` now publishes the `.faa` too, but the GFF3 route stays:
the coordinates, strand and phase Pyrodigal called, translated against the same
window, reproduce the proteins exactly.

### Clustering engine

| engine | when | why |
|---|---|---|
| `blastp` | the default here | already installed natively; 363 proteins is nothing |
| `mmseqs` | the NCBI sweep | all-vs-all blastp does not scale past a few thousand |

`--engine auto` uses `mmseqs` if it is on `PATH` and `blastp` otherwise. Both
apply the same cuts — 90% identity, 80% coverage **of both sequences** — which
is the intent of Acman's Roary `-i 0.9 -s`.

MMseqs2 is not installed on this laptop or on Gadi. The official universal
macOS build (native arm64, no Homebrew, no conda) is:

```bash
mkdir -p ~/opt/mmseqs2 && cd ~/opt/mmseqs2
curl -LO https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-osx-universal.tar.gz
tar xzf mmseqs-osx-universal.tar.gz
export PATH="$HOME/opt/mmseqs2/mmseqs/bin:$PATH"
```

### What it produced

7 loci, 363 CDS → **258 families**, of which **57 span more than one locus**.
The 140 `hypothetical protein` CDS that shared one colour split into **107
families**.

The most widely shared families are the ones that should be, which is the
check that the clustering is sane:

| family | loci | aa | product |
|---|---|---|---|
| `Feff4f214` | 5 | 79 | RelB |
| `F7b88ed79` | 5 | 95 | RelE |
| `F9c5dc466` | 4 | 446 | ISEsa1 transposase |
| `F97780e40` | 4 | 307 | Cation-efflux pump (czc) |
| `F5dbb1fdc` | 4 | 412 | ISEsa2 transposase |
| `F2bc0d751` | 4 | 190 | Cytochrome b561 (cytb) |
| **`F61fefa1c`** | **4** | **70** | **hypothetical protein** |
| **`F0e2e43c1`** | **4** | **55** | **hypothetical protein** |

The last two are a result. Two unnamed families sit in 4 of the 7 loci —
as conserved as czc and cytb, and invisible until now because no database
could name them. They are candidates for the conserved core that Stage 2 has
to delimit.

Family ids are a hash of the representative sequence, not a cluster index, so
a family keeps its identity and its colour as loci are added.

---

# STEP 12 — Draw the figure

This one runs **on the laptop**, not on Gadi: it is the iterate-and-look step,
it needs no database, and matplotlib is already installed at
`/opt/anaconda3/bin/python3.12`. Copy the three tables down first.

```bash
python3 scripts/colour_blocks.py \
  --blocks   run_bakta/blocks \
  --context  run_bakta/context.tsv \
  --metadata run_bakta/record_metadata.tsv \
  --families results/gene_families.tsv \
  --outdir   figures
```

Drop `--families` and the CDS track goes back to flat grey.

### STEP 13 — Group isolates by shared history

```bash
python3 scripts/structural_groups.py \
  --blocks run_all/blocks \
  --out    results/structural_groups.tsv \
  --out-pairs results/structural_pairs.tsv
```

Jaccard similarity on which **vehicle + core context genes** each locus
carries -- `ISEsa1, ISEsa2, cytb, czc, relB, relE, hok_gef` by default (the
same five names `draw_schematic()` in `colour_blocks.py` already calls "the
conserved core"), **not** every IS and cargo call in the ±20 kb window.

The first version of this script scored the whole window and found almost
nothing: even the six COMPLETE loci, known to share one vehicle at one
offset, scored 0.43-0.57 against each other. The variable upstream and the
`tra` tail are -- by this project's own finding -- "a different history in
every locus", so including them in a score meant to find SHARED history works
against it: one locus-private IS anywhere in 40 kb costs real Jaccard when
the feature set is this small. Restricting to the vehicle and core genes that
are actually meant to be shared fixed it.

A fragment does not count as present. A 64 bp or 275 bp stub is evidence an
event happened, not evidence the locus currently carries a working copy, and
counting it as "present" would blur exactly the COMPLETE/DECAYED distinction
this grouping exists to sharpen.

**Result at the default 0.90 threshold, 14 loci:**

| group | loci | shared |
|---|---|---|
| 6 | CP012599, CP018961, CP020841, CP021960, CP028197, CP028537 | ISEsa1, ISEsa2, cytb, czc, hok_gef, relB, relE |
| 3 | CP023893, CP026975, LT991956 | ISEsa1, hok_gef, relB, relE |
| 2 | CP006642, CP017220 | relB, relE |
| singleton | CP010512 | ISEsa1, cytb, czc, hok_gef, relB, relE |
| singleton | CP026216 | ISEsa1, cytb, czc, relB, relE |
| singleton | MF344582 | ISEsa2, cytb, czc, hok_gef, relB, relE |

Every multi-member group is **pairwise 1.00 internally** -- not merely above
threshold, identical. The three singletons each carry a genuinely different
subset (CP010512 has no ISEsa2; CP026216 has no hok_gef; MF344582's only
complete vehicle copy is ISEsa2, since its ISEsa1 is the 275 bp fragment) and
correctly fall outside every group rather than being forced into the nearest
one.

Grouping is **union-find at the threshold**: two loci are linked iff their
core Jaccard clears `--min-similarity`, and groups are connected components.
This can chain (A-B and B-C both clear it without A-C clearing it) -- the
output's `chained` column flags any group where that happened. It did not
happen here; every multi-member group's `weakest_internal_pair` is 1.00.

**What this is not.** It is the fast triage CLAUDE.md Stage 6 calls for
before the sequence-level pass, not a replacement for the decisive test in
Stage 5 (all-vs-all `blastn -task dc-megablast`, threshold-stepped). Two loci
could share an element catalogue by convergence rather than by one insertion
event -- this script cannot tell those apart, because it never looks at the
DNA between the elements, only at which named elements are present.

Colour is **content-addressed on the shared element set**, not on membership,
so a group keeps its colour across runs as its architecture signature
repeats -- the convention `colour_blocks.py` already uses for IS-family and
gene-family colour. Singletons share one neutral grey rather than each
getting a near-identical hue nobody could tell apart.

### Showing the groups on the figure

```bash
python3 scripts/colour_blocks.py ... \
  --groups results/structural_groups.tsv --outdir figures
```

Adds a coloured swatch beside each locus label and a dedicated row under the
figure listing each group's members and shared elements. Loci with no
multi-member group (the three singletons above) get no swatch.

### Styles

`--style simple` is the default: **one row per locus**, IS and cargo on the
same line, colour encoding **role and orientation only**.

| | colour |
|---|---|
| IS, forward | mustard `#d4a017` |
| IS, reverse | yellow `#f5dd5d` |
| relB | green `#00b050` |
| relE | red `#e8112d` |
| cytb (cytochrome b561) | purple `#7030a0` |
| czc (Co-Zn-Cd efflux) | blue `#2e75b6` |
| ars (arsenic resistance) | teal `#00b0a0` |
| tra (conjugation) | dark grey `#4a4a4a` |
| sil, intI, qnr, hok_gef | light grey `#c9c9c9` — labelled, not singled out |
| CDS | near-white, thin underlay |

Which IS it is comes from the **label**, not the colour, so the two IS tones
are free to carry orientation instead. That is the whole trade: you lose
family at a glance and gain strand at a glance, and the figure answers one
question rather than twenty-four.

CDS stay as a thin pale underlay rather than disappearing, so the backbone is
still visible and an IS can be seen to contain CDS rather than replace them.

`--style tracks` keeps the detailed version — IS, CDS and cargo on separate
lines, IS coloured by family with the element as a shade, CDS coloured by gene
family. Use it when the question is *which* element or *which* family; use
`simple` when the question is shape.

Cargo blocks are labelled on the figure under the same deduplicating rule as
IS, so `czc`, `cytb` and `ars` can be read directly rather than matched back
to a legend swatch.

### Results at 14 loci (`REPLICON=all`)

The full run adds the seven chromosomal loci. The seven plasmid rows in
`context.tsv` are **byte-identical** to the plasmid-only run, so nothing about
adding them perturbed the existing calls.

| | COMPLETE | DECAYED | HALF-CORE | NO CORE |
|---|---|---|---|---|
| plasmid (7) | 3 | 1 | 2 | 1 |
| chromosome (7) | 3 | 0 | 3 | 1 |

Consensus offsets, now over 14:

```
element   loci   offset   spread
ISEsa1      11     -5494   -12244 .. -3991
ISEsa2       7     -3898    -6675 .. -3897
cytb         9     -1960    -4738 .. -1960
czc          9     -1184   INVARIANT
relB        14        -2        -2 .. 2
relE        14      +241       241 .. 5290
hok_gef     11      +596       595 .. 1697
```

**czc is at exactly −1184 in all nine loci that carry it** — four plasmids and
five chromosomes, across five genera. Doubling the dataset did not loosen it.
That is the strongest candidate yet for the left edge of the conserved core
that Stage 2 has to delimit.

#### The one-sided asymmetry is scale-dependent

Slide 11 states that every insertion is upstream of the core and only KPN1482
shows anything downstream. Measured across all 14 at ±20 kb:

| | loci |
|---|---|
| IS **within +3 kb** of the bait | **1** of 14 (CP006642, ISSen3 at +762) |
| IS **beyond +3 kb** | 5 of 14 |
| upstream only | 9 of 14 |

So the claim holds at the scale the deck drew — its panels are ~7–8 kb wide —
and breaks at ours. The honest statement is **the region immediately
downstream of relE is near-inviolate; the region beyond +3 kb is not**.

CP020841 and CP021960, two independent *K. pneumoniae* chromosomes, both carry
ISEcp1 at **exactly +3436**. An identical downstream offset in two records is
a shared event, not two insertions.

### Window size: keep ±20 kb, zoom the figure

A ±5 kb window is tempting — the conserved core spans only −1960 to +596 — but
it would destroy the vehicle call. Measured against the 14 loci:

| ISEsa1 at ±5 kb | loci |
|---|---|
| lost entirely | **4** (CP023893, CP026975, LT991956, and CP018961, *both* pieces) |
| clipped at the window edge | 6 |
| intact | **2** of 11 |

CP018961 would lose ISEsa1 *and* ISEsa2 and report NO CORE — inverting the
Stage 1 reference result. ISEsa1's own median offset is −5494, already outside
a ±5 kb window.

The extraction window and the figure's x-range are different things, so keep
the window at ±20 kb and zoom the **figure**:

```bash
python3 scripts/colour_blocks.py ... --xlim=-8000,3000 --outdir figures_zoom
```

`--xlim` changes nothing but the view: no re-extraction, no re-analysis, no
call changes. Use the `=` form — argparse reads a bare `-8000,3000` as a flag.

Two things it does so the zoom cannot mislead:

- **features outside the range are clipped, not parked at the edge.** Without
  that, an element at −12 kb is drawn hard against the axis in a figure zoomed
  to −8 kb, on top of the locus labels, looking like data that is off-screen.
- **the title counts what is hidden** — `74 IS/cargo outside`. A figure that
  silently drops a third of its calls looks like a finding.

Labels anchor to the display range rather than the extraction flank. They used
to anchor to the flank, which put the locus names 20 kb off-screen in a zoomed
figure and made `bbox_inches='tight'` grow the canvas to include them: a
6161 px image of an 11 kb window.

### The figure states its own scope

The title carries the replicon filter — `[plasmids only — 7 loci]`. A figure
showing 7 of the 14 records that had a bait hit, without saying which 7, would
overstate what it covers. If every locus drawn shares one `replicon_type`,
that is a **filter, not a finding**, and the title says so.

The chromosomal loci on slides 2 and 4 of the 2018 deck are not missing from
the data. They were found and then dropped by `REPLICON=plasmid`, which is
`env/gadi_run.sh`'s default:

```
replicon_type    chromosome    7   dropped
replicon_type    plasmid       8   kept
```

To include them:

```bash
qsub -v BAKTA=true,BAKTA_DB=/g/data/rg42/bakta_db/db-light,\
REPLICON=all,OUT=/scratch/rg42/$USER/run_all env/gadi_run.sh
```

### The idealised-architecture panel

Below the loci, separated by a rule, is a consensus schematic. **It is computed
from the data, not asserted.** Each block sits at the **median offset** of that
element across the loci carrying it, complete copies only — a fragment's offset
is where what survived begins, not where the element was, so averaging the two
would put the consensus somewhere neither has ever been.

```
  element   loci   offset   spread
  ISEsa1       5    -10554   -12244 .. -5494
  ISEsa2       4     -3898   -6675 .. -3898
  cytb         4     -1960   -4738 .. -1960
  czc          4     -1184   INVARIANT
  relB         7        -2   -2 .. 2
  relE         7      +241   241 .. 5290
  hok_gef      6      +600   596 .. 1697
```

An element at an identical offset in every locus carrying it is marked
`(same in all N)` and set in bold. Only **czc** earns that here — −1184 in all
four loci that have it — and relB is within 4 bp of invariant. That is the
difference between a schematic and a claim: czc is drawn at −1184 because the
data put it there in every case, not because a slide says so.

The panel also shows what is *not* fixed: a dashed band upstream for the
variable IS region, and a `tra` band spanning its full observed range
(+9246 to +19753) rather than a single position, because it has one.

Turn it off with `--no-schematic`.

One thing to watch: `sil, intI, qnr, hok_gef` and `CDS` are both greys and are close at
small sizes. They are distinguishable by height (cargo full, CDS thin), but if
that matters for a given panel, change `SIMPLE['cargo']` in the script.

### Labelling — one rule, not three special cases

Every IS call gets a name on the figure. The rule used to be
`name.startswith('ISEsa') or fragment or nested`, which had this project's
vehicle hardcoded into it: run the pipeline with a different bait and the
elements that mattered would have gone unlabelled with no warning — the same
failure `--vehicle_ends` was added to `classify_context.py` to prevent.

| flag | effect |
|---|---|
| `--label-is distinct` | **default in `simple`.** One label per distinct element; a fragment counts separately from a complete copy, so CP018961's 1540 bp ISEsa1 and its 64 bp remnant are both named. |
| `--label-is all` | default in `tracks`. Every IS named, 45 of 45 in this set. |
| `--label-is key` | only the vehicle ends, fragments, and nested elements |
| `--label-is none` | no labels |
| `--vehicle-ends ISEsa1,ISEsa2` | which elements are **bold**. Must match the `ENDS` the run used, or the figure emphasises the wrong ones. |
| `--no-label-names relB,relE` | never labelled, because their colour is already unique. Pass `''` to label everything. |

`relB` and `relE` are green and red and sit at the same place in every row, so
a label adds nothing — and it cost something. Labelling them bold forced two
extra stack levels at x=0 on all 14 rows, and row pitch is derived from the
deepest stack, so two redundant labels were setting the height of the whole
figure. Dropping them took it from 4508 px to 3055 px.

They keep their labels in the idealised-architecture panel, where the text
carries the **offset** (`relB -2`, `relE +241`) rather than the identity.

Labels that would overprint are stacked upward, each with a leader line, and
**the row height is computed from the deepest stack** rather than fixed — so a
locus with thirteen IS calls (CP006642, three levels) does not push its labels
into the row above. Length is appended only where it is the evidence: on
fragments.

Verify coverage with:

```bash
python3 - <<'EOF'
import csv, glob, sys; sys.path.insert(0, 'scripts')
from colour_blocks import label_plan
for f in sorted(glob.glob('run_bakta/blocks/*.blocks.tsv')):
    rows = list(csv.DictReader(open(f), delimiter='\t'))
    n = sum(1 for r in rows if r['feature_type'] == 'IS')
    p = label_plan(rows, 'all', {'ISEsa1', 'ISEsa2'})
    print(f, n, len(p), max(p.values(), default=0) + 1)
EOF
```

Cargo is deliberately **not** labelled on the figure — twenty `tra` blocks
would bury everything else — and is read from the legend instead.

### Palettes

`--palette` picks the colour scheme. Both deck palettes were read out of the
`.pptx` shape XML, not sampled from a screenshot, so the hex values are the
decks' own.

| `--palette` | source | relB | relE | ISEsa1 | ISEsa2 | core (cytb, czc) |
|---|---|---|---|---|---|---|
| **`deck`** (default) | `relBE_schematic.pptx` slide 11 | `#028090` teal | `#d1495b` red | `#5ba4a8` teal | `#b57f45` brown | `#d8e0e0` pale |
| `lab2018` | Sept 2018 lab meeting, slide 2 | `#00b050` green | `#ff0000` red | `#00b0f0` blue | `#ffc000` amber | white |
| `default` | this script's own | blue | cyan | red | orange | green shades |

`deck` is the default because slide 11 is the nine-loci figure this pipeline
reproduces, and CLAUDE.md names that deck as the source.

IS families the decks never drew are given tints **constrained to the deck's
warm hue range**, so the variable upstream still reads as one kind of thing at
a glance — which is the decks' point — while staying separable close up.

One thing carried over faithfully that you may not want: both decks paint
`Cytb/b6` and `CO-Zn-Cad` the **same** colour, because they are presented as a
single conserved core rather than as two genes. In this figure they are
adjacent blocks and therefore hard to tell apart. Give them distinct values in
`PALETTES` if you would rather read them separately.

`run_bakta/` is the Bakta-annotated run pulled down from Gadi:

```bash
scp -r gadi:/scratch/rg42/$USER/run_bakta/results/{blocks,context.tsv,record_metadata.tsv} run_bakta/
```

Point `--blocks` at `results/blocks` instead to redraw from the GenBank
annotation — the `annot_source` column in the blocks table says which you have.

Three files land in `figures/`:

| file | what it shows |
|---|---|
| `colour_blocks.png` / `.svg` | one row per locus, x = distance from the bait. IS track above each line, curated cargo below, CDS as grey context. |
| `is_matrix.png` | loci x IS element presence, columns grouped by IS family |
| `is_matrix.tsv` | the same matrix as a table — this is the input to the Jaccard/PERMANOVA grouping |

The `.svg` is there so panel labels and fonts can be edited in Illustrator or
Inkscape without re-running anything.

## Reading it

**Vertical position is the track, not the strand.** IS above each line, CDS on
it, cargo below. This was misread as strand the first time someone looked at
the figure, so strand is now drawn explicitly and in only one place: **the
direction each block points**.

In CP018961, for example, `relB` and `relE` are drawn *below* the line because
they are cargo, and they point *right* because they are on the plus strand:

```
cargo   relB   strand=+   offset=-2
cargo   relE   strand=+   offset=241
```

The `minus` in the locus id `CP018961__CP018961__1__minus` refers to the
orientation of relBE in the **original record**, before the window was flipped.
Every window is reoriented so relB is on `+` (CLAUDE.md, "Reorient first"), so
within the figure's coordinates relB is always on the plus strand — if it were
not, every offset in every table would be meaningless. Check the FASTA header:

```
>CP018961__CP018961__1__minus orig=CP018961:237696-278222 strand=minus ...
```

`strand=minus` there is the original; the sequence that follows is the flipped
one.

**IS colour is the family; shade is the element inside it.** Twenty-four
elements given twenty-four unrelated colours is a lookup table, not a picture.
ISfinder names already carry the family (`ISEsa1_ISNCY_IS1202` → ISNCY), so
twelve hues now do the work of twenty-four, and related elements look related.

Two things that per-element colour hid, and this does not:

- **ISEsa1 is ISNCY; ISEsa2 is IS110.** The two ends of the candidate vehicle
  are *different families*. The composite-transposon rule in CLAUDE.md Stage 3
  looks for two IS of the **same** family flanking the cargo — it does not
  describe this structure, and should not be expected to fire on it.
- **ISEsa1 and ISKpn21 are the same family and group** (ISNCY/IS1202). They now
  share a hue, which is exactly the pair that makes a naive same-family rule
  fire on CP018961 — where ISKpn21 is nested *inside* ISEsa1 and flanks nothing.
  This is the false positive recorded as finding 4 in the project plan, and the
  figure now shows the reason for it rather than hiding it in a name.

**The vertical dotted line at −3898** is the claim the figure exists to test.
ISEsa2 begins at exactly that offset in CP028197, CP028537 and MF344582 —
three records, three different genera. An offset that recurs to the base pair
across independent loci is the signature of **one ancestral vehicle**, not of
three separate insertions that happened to land nearby.

**One element can appear twice without being two insertions.** In CP018961:

```
ISEsa1   ISNCY   -10554  1540 bp  complete
ISKpn21  ISNCY    -9017  2278 bp  complete   nested_in=ISEsa1
ISEsa1   ISNCY    -6739    64 bp  fragment
ISEsa2   IS110    -6675  1492 bp  complete
```

That is **one** ISEsa1 copy — 1540 + 64 = 1604 bp, its full length — split in
two by ISKpn21 landing inside it. The figure draws a bracket over the pieces
labelled `one ISEsa1, 1604 bp in 2 pieces` so it cannot be misread as two
independent insertions.

The pale block between them is **ISKpn21**, not a second colour of ISEsa1.
ISKpn21 is the same family as ISEsa1 (ISNCY/IS1202) and therefore the same
hue, one shade lighter. It is labelled now: any element nested inside another
gets a label, because that nesting *is* the structural event in its locus, and
leaving it blank made the row read as "ISEsa1 drawn twice in two colours".

This pair is also the composite-transposon false positive from the project
plan — same family, same group, but one is inside the other and flanks
nothing.

**Hatched blocks are fragments** (<50% of the element present). The 64 bp
ISEsa1 remnant in CP018961 is the smallest and the most important: it sits at
−6739 and ends at −6675, which is precisely where that locus's ISEsa2 starts.
This is the direct analogue of the ~75 bp ISAba125 remnant that identified
Tn125 as the ancestral vehicle for blaNDM. It is only in the figure because
`is_search.py` runs with **no length cutoff** — a 50 bp or 100 bp minimum,
which most IS pipelines apply by default, would delete it silently.

**A dashed block outline** means the element runs off the edge of the ±20 kb
window and is not actually a fragment. Do not read those as decay; the
`clip_cause` column in `results/is/*.tsv` says which is which.

## What Bakta changed, measured

Re-annotating with Bakta `light` and redrawing (gate passed: `context.tsv`
byte-identical, so nothing leaked into the BLAST layers):

| locus | GenBank CDS / hypo | Bakta CDS / hypo | |
|---|---|---|---|
| CP006642 | 32 / 7 (22%) | 30 / 13 (43%) | **worse** |
| CP018961 | 39 / 19 (49%) | 38 / 19 (50%) | flat |
| CP023893 | 60 / 27 (45%) | 48 / 11 (22%) | much better |
| CP028197 | 51 / 33 (65%) | 47 / 22 (46%) | better |
| CP028537 | 50 / 32 (64%) | 49 / 20 (40%) | better |
| LT991956 | **0 / 0** | **57 / 22 (38%)** | from nothing |
| MF344582 | 46 / 35 (76%) | 45 / 20 (44%) | much better |
| **total** | 278 / 153 (**55%**) | 314 / 127 (**40%**) | |

Two things worth keeping in view.

**LT991956 went from contributing nothing to contributing 57 CDS.** It is the
one record with no CDS features at all, so before this it was a blank line in
the figure. That alone justified the run.

**CP006642 got worse** — 7 unnamed became 13. GenBank's curated annotation of
that record names more products than Bakta's `light` database can. This is the
first concrete evidence that `full` would buy something, and it is confined to
one locus; it does not change any coordinate, any IS call or any cargo hit.

**CDS colour is the gene family**, from STEP 11b. There is deliberately no
legend for it: 57 colours cannot be told apart from a key, and the figure does
not ask you to identify a family. It asks whether the **same colour appears at
the same offset in another row** — which is the whole Acman Figure 4 device.
Grey means the family occurs in only one locus (`--min-family-loci` controls
the cut; 1 would give every singleton its own colour and turn the track into
confetti).

## What the figure does NOT yet show

Absent, and absent on purpose: the splitting tree that runs down the left
of Acman Figure 4. That needs Stage 5's all-vs-all threshold stepping. With
seven loci it would be a tree of seven leaves and would imply far more support
than seven plasmids can carry.

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
