# TAflank — implementation description

How the pipeline is built, what each step does, how hits are filtered, how
flanking sequence is retrieved, and what comes out.

Companion documents, kept deliberately separate:

| file | holds |
|---|---|
| `CLAUDE.md` | the project plan and conventions |
| `PIPELINE.md` | the **scientific rationale** — why each filter exists (F1–F6), the MGE layers, the decision tree |
| **`PIPELINE_DESCRIPTION.md`** (this file) | the **implementation** — what actually runs, with what parameters, producing what files |

---

## 1. What it does

Takes a **bait** sequence, finds every copy of it in a set of genomes, works out
whether each copy sits on a plasmid or a chromosome, cuts the neighbourhood out
around it, and characterises the mobile genetic elements in that neighbourhood.

The bait is a **parameter**. Nothing in the pipeline hardcodes relBE — point
`--bait` at any gene and the whole chain works.

```
bait.fa ──► local BLAST ──► Tier A ──► fetch metadata ──► Tier B ──► ±20 kb windows
 (any gene)   (vs a DB)    (sequence)    (batched)        (record)        │
                                                                          ├─► IS elements (ISfinder)
                                                                          ├─► replicon type (PlasmidFinder)
                                                                          └─► collated tables
```

## 1b. Importing an NCBI web BLAST result

A web BLAST sees a far larger database than any local one — `nt` is
3,697 Gbp against `public16`'s 0.047 Gbp, about 79,000x — so a search there
legitimately returns thousands of hits where a curated local set returns
fifteen. Those results can be fed straight in:

```bash
# On the NCBI results page: Download -> Hit Table (text)  [or csv]
./tools/nf run main.nf -profile local \
    --hit_table ncbi_hits.txt --genome_dir ''
```

`BUILD_DB` and `BAIT_SEARCH` are skipped; `IMPORT_BLAST` takes their place and
applies **the same Tier A thresholds**, so an imported search and a local one
are directly comparable.

Two consequences follow from what a hit table does and does not carry:

| column | local search | hit table |
|---|---|---|
| the 12 `outfmt 6` fields | yes | **yes — a hit table *is* outfmt 6** |
| `qlen` | yes | taken from the bait FASTA |
| `slen`, `stitle` | yes | **absent** |

- **`meta_source` is forced to `ncbi`.** The `blast` mode parses `stitle`,
  which a hit table does not have, so metadata comes from `esummary` instead.
- **Windows are fetched from NCBI**, by accession and range, since there is no
  local copy of the record. Only the needed span is requested, except where a
  circular window wraps the origin — that needs the whole record for the
  arithmetic. Responses are cached under `data/efetch_cache/`.
  PlasmidFinder fetches the **whole** replicon, because Inc *rep* genes are
  rarely within 20 kb of the bait.

Verified: windows fetched this way are **byte-identical** to locally extracted
ones, and report true replicon coordinates
(`orig=CP018961.1:237696-278222`), not coordinates within the fetched slice.
Vehicle classification agrees locus for locus.

## 2. Quick start

```bash
./tools/nf run main.nf -profile local                      # plasmids (default)
./tools/nf run main.nf -profile local --replicon_type chromosome
./tools/nf run main.nf -profile local --replicon_type ''   # everything
./tools/nf run main.nf -profile local --bait anchors/other_gene.fa --subregions ''
```

On a cluster against a pre-formatted database, with no network at all:

```bash
./tools/nf run main.nf -profile gadi \
    --blastdb /scratch/$PROJECT/$USER/blastdb/nt_prok \
    --genome_dir '' --meta_source blast --taxonomy
```

Use `./tools/nf`, not bare `nextflow` — see §8.

---

## 3. Stages

### BUILD_DB — `scripts/build_blastdb.sh`

Concatenates a directory of FASTAs and builds a nucleotide BLAST database with
`-parse_seqids -hash_index`.

`-parse_seqids` is not cosmetic. It is what makes **§6 sequence retrieval by
accession and coordinate range** possible, and therefore what lets the pipeline
scale past the point where you can keep every genome as a loose file.

Emits a `.dbmeta.tsv` sidecar — source path, sha256, sequence count, total bp,
BLAST version, build time. This exists because two earlier databases
(`blastdb/all16`, `blastdb/tadb3_loci`) recorded their source as `/tmp` files
that no longer exist, making them untraceable and, in one case, silently stale.

### BAIT_SEARCH — `scripts/bait_search.py`  *(Tier A)*

`blastn` of the bait against the database. **`-task blastn`** (word size 11),
deliberately not megablast, which misses diverged copies at 85–92% identity.

`--subregions relB:1-240,relE:240-527` names parts of the bait so that a
half-length hit can be called an orphan antitoxin or an orphan toxin rather
than lumped as "partial". Omit it for a bait with no meaningful internal parts.

Emits `sacc`, not `sseqid`: with `-parse_seqids`, `sseqid` renders as
`gb|CP018961.1|` and breaks every downstream join on accession.

### FETCH_METADATA — `scripts/fetch_metadata.py`  *(Tier B input)*

Three sources, chosen with `--meta_source`:

| mode | source | network |
|---|---|---|
| `ncbi` | batched `esummary`, cached per accession | yes, once |
| `cache` | the on-disk JSON cache | no |
| `blast` | the `stitle` and `slen` columns of the BLAST output | **no** |

`blast` mode is what makes the search runnable on a compute node with no
outbound access, and it is the only workable mode for a large sweep, where one
`esummary` per hit would be hopeless. Agreement with `esummary` on the
16-record set is **14/14**.

**Replicon classification** uses the structured `genome` field where present.
That field is not always populated — CP021960 is a 5.6 Mb circular
*K. pneumoniae* "complete genome" with `genome: None` — so it falls back to the
submitter's title and records which route was taken in
`replicon_type_evidence` (`esummary_genome` | `title` | `none`).

Before any chromosome rule there is a **draft-assembly guard**: titles
containing *whole genome shotgun*, *WGS*, *scaffold*, *contig*, *unplaced* or
*unlocalized* return `unknown`. `nt_prok` is full of these and they are
genuinely neither a plasmid nor a chromosome. The guard runs first because such
titles often still contain the word "genome".

There is deliberately **no size heuristic**: megaplasmids exceed 1 Mb and would
be misassigned as chromosomes.

`replicon_type_evidence` records which rule fired — `esummary_genome`,
`title_plasmid`, `title_chromosome`, `title_complete_genome`,
`title_complete_sequence`, `title_draft_assembly` or `none` — so a structured
call is always distinguishable from an inferred one.

Records NCBI cannot resolve get `replicon_type=unknown`, never forced into a
bucket. PlasmidFinder gives independent evidence for those — a replicon hit
says "plasmid" even when the metadata is silent.

### SELECT_LOCI — `scripts/select_loci.py`  *(Tier B)*

Joins hits to metadata, computes flank availability, applies the filters, and
marks each hit `selected` yes/no with a `deselect_reason`. See §4 and §5.

### EXTRACT_LOCI — `scripts/extract_loci.py` → `scripts/relbe_window.py`

Cuts the ±`flank` window and reorients so the bait is on the `+` strand. See §6.

### IS_SEARCH — `scripts/is_search.py`  *(PIPELINE.md layer 2 + F5)*

`blastn` of each window against the ISfinder mirror, **with no length cutoff**.

> `--min-length` defaults to 0 and the script **asserts** it is 0. The whole
> ancestral-vehicle test depends on a ~64 bp ISEsa1 remnant, mirroring the
> ~75 bp ISAba125 remnant that identified Tn125 for *bla*NDM. A length cutoff
> anywhere in this layer destroys the key evidence.

**Tiling** keeps the best-scoring HSPs, allowing `TOLERANCE = 30 bp` of overlap.
Never zero: PIPELINE.md §4 records a real incident where a 3 bp overlap at zero
tolerance discarded a 1,540 bp ISEsa1 alignment, making a near-complete element
look like a bare stub and inverting the biological reading.

**Coverage class** by subject coverage — `complete` ≥90%, `partial` 50–90%,
`fragment` <50%. Fragments are kept; they are the evidence.

**Window-edge clipping** (`clipped_at_window_edge`, `clip_cause`,
`clipped_missing_bp`). A low `query_cov` has two unrelated causes, and the
decisive test only cares about one of them:

- the element is a **genuine decayed remnant** sitting in open DNA — evidence;
- the element is **whole in the real replicon** and was merely sliced by the
  edge of the ±20 kb window *we chose* — an artefact of construction.

Both arrive from `blastn` as `is_fragment=yes`. Counted together they would put
phantom remnants at exactly ±20000 in the ISEsa1 offset distribution — the one
plot the ancestral-vehicle conclusion rests on.

An alignment touching an edge is only called clipped if part of the subject is
left unaccounted for **on that same side**, i.e. the element demonstrably runs
on past where we stopped looking. Which subject end that is depends on strand:
on `-` the subject runs backwards relative to the window.

`clip_cause` keeps the two causes apart, and the distinction is not cosmetic:

| `clip_cause` | meaning | how to read it |
|---|---|---|
| `window` | the ±20 kb cut ended there while the replicon continued | **discount as evidence** — widen the window to resolve |
| `replicon` | the replicon or contig genuinely ends there (header says `truncated_up/dn=yes`) | **real truncation** — a signal, per CLAUDE.md |
| `mixed` | one edge of each | both |

Nothing is filtered. `is_fragment` keeps its plain meaning (`query_cov < 50`)
and the new columns say whether that number can be read as decay.
Downstream, `n_is_fragment_real` (context, summary) and
`n_is_fragment_clipped` (summary) separate the counts, and
`vehicle_edge_clipped` names any ISEsa1/ISEsa2 end whose true coverage is
**unknown** rather than low — such an end must never drive a `DECAYED` call.

`cargo_search.py` carries the same three columns for the same reason: a gene
lying across the boundary reports low `cov` because we stopped looking, not
because the gene is broken.

**Nesting** (`nested_in` / `interrupts`) requires *both*:
1. element B lies within the window span of element A, and
2. A's own HSPs are **near-contiguous in subject coordinates** and sum to no
   more than one copy length.

Condition 2 is essential. Without it, several independent full-length copies of
one element are misread as a single split copy and everything between the
outermost two is falsely flagged as nested — which fired 8 false events on the
IS26/IS15DI array in CP006642.

### PLASMIDFINDER — `scripts/plasmidfinder.py`  *(layer 5)*

Replicon (Inc) typing by plain `blastn` against the CGE PlasmidFinder FASTA.
No conda, no CGE wrapper, no Docker — the database is plain FASTA and the
native BLAST types it directly.

Thresholds are Acman et al. 2022: **coverage ≥80%, identity ≥95%**, with
overlapping hits on one replicon resolved in favour of the higher identity.

Two things that are easy to get wrong and are handled explicitly:

- Runs on the **whole parent replicon, not the window**. Inc *rep* genes are
  rarely within 20 kb of the bait, so a window-only run returns empty for
  nearly everything.
- The replicon gene is the **query**, the replicon the **subject**, so coverage
  is `qcovs` — coverage of the *rep gene*. Getting this backwards is the same
  `qlen`-vs-`slen` trap that once misclassified all nine loci.

> **Taxon scope.** PlasmidFinder ships databases for Enterobacteriales and
> Gram-positives only. Records outside that scope must report
> `NA_no_db_for_taxon`, **never `none`** — absence of a database is not absence
> of a replicon. Relevant as soon as the TADB3 set (Cyanobacteria,
> Actinobacteria, *Helicobacter*, *Caulobacter*) is added.

### CARGO_SEARCH — `scripts/cargo_search.py`  *(Tier C input)*

`blastn` of each window against `assets/cargo_genes.fa` — 333 sequences across
12 families (`czc`, `cytb`, `hok_gef`, `tra`, `qnr`, `intI`, `mer`, `ars`,
`sil`, `ter`, `relB`, `relE`), built by `scripts/build_cargo_set.py` from this
project's own GenBank records, so every entry traces to an accession and
coordinates.

Deliberately **not** an annotation pass. Classifying context needs only the
families the biology turns on, and 333 sequences answer that in milliseconds.
Same 30 bp tiling tolerance as the IS layer.

### CLASSIFY_CONTEXT — `scripts/classify_context.py`  *(Tier C)*

Combines the IS and cargo layers into `results/context.tsv`: a vehicle class
(COMPLETE / DECAYED / HALF-CORE / NO CORE) plus a cargo profile and the offsets
of `czc` and `hok_gef`.

`vehicle_interrupted` is reported separately from the class, because
"all sequence present" and "intact" are different things — CP018961 is COMPLETE
at 96% ISEsa1 coverage yet split in two by a 2,278 bp ISKpn21 insertion.

> Coverage is judged on the **best HSP per element**, not per HSP. An element
> split by a later insertion contributes a large piece and a small one; judging
> on the small one would call a 96%-covered ISEsa1 "decayed" purely because its
> 64 bp tail is also reported.

Reproduces the hand classification on all four overlapping loci.

### BAKTA + GENE_BLOCKS — `scripts/gene_blocks.py`

Bakta and ISfinder are **complementary and independent**, so they run in
parallel and are merged afterwards:

```
            ┌─► IS_SEARCH  (ISfinder)   → names the ELEMENT   ISEsa1
locus ──────┤
            └─► BAKTA      (annotation) → names the FAMILY    "ISNCY family transposase"
                      │
                      └──► GENE_BLOCKS  ← element name wins on >=50% overlap
```

**The two-tier colour key is not cosmetic.** Annotation names IS elements by
family only, so colouring from it alone would merge the vehicle with unrelated
transposases of the same family. Measured on CP028197:

| annotation says | must say |
|---|---|
| "ISNCY family transposase" | **ISEsa1** |
| "IS110 family transposase" | **ISEsa2** |
| "IS91 family transposase" | **ISVsa3** |
| "IS6 family transposase IS15DIV" | **IS26** |

Annotation source falls back **GFF → parent GenBank → none**, recorded per
locus in `annot_source`. The GenBank fallback is not a shortcut: it lets the
merge be built and validated before Bakta exists, since Bakta cannot run on
macOS at all. Swapping Bakta in changes the source, not the logic.

The gap stays visible — `LT991956` yields **0 CDS** from GenBank, so its gene
track is empty and the table says so. That locus is the concrete argument for
running Bakta on the cluster.

`--run_bakta` defaults to false; `GENE_BLOCKS` joins on `remainder: true`, so
every locus gets a track whether or not Bakta ran.

### COLLATE — `scripts/collate.py`

Gathers per-locus and per-record outputs into the deliverables, and attaches
`lab_alias` from `assets/replicon_aliases.tsv` so tables carry **your** names as
well as NCBI's — the deck calls LT991956 `JIE4345` and MF344582 `p525011-HI2`,
where NCBI calls them `pC309-p2` and `pNY2385-HI2`.

Every locus in `loci.tsv` appears in `summary.tsv` even if a layer produced
nothing for it, so a silent drop shows up as a visible gap rather than a
quietly smaller number.

---

## 4. On what basis hits are filtered

Three tiers. Numbers are from the default plasmid run on the 16-record set.

```
blastn bait vs DB                                            15 HSPs
      │
TIER A  sequence only — free, from the BLAST output alone
      │    FULL        pident ≥95 AND qcov ≥95                    14
      │    ORPHAN_*    half hit, one named subregion ≥80% covered  0
      │    PARTIAL     below qcov, no single subregion dominant    0
      │    PARALOGUE   pident <85 → RelE/ParE superfamily          1
      ▼
TIER B  record metadata — batched esummary
      │    replicon_type = plasmid        → dropped 7 chromosomes
      │    tier_a        = FULL           → dropped 1 paralogue
      │    flank availability             → min_flank=0, dropped none
      │    is_complete                    → not required by default
      ▼
                                                              7 selected
TIER C  context — needs the window (not yet implemented)
           does the ISEsa1/ISEsa2 vehicle co-occur?
           → COMPLETE / DECAYED / HALF-CORE / NO CORE
```

**Why Tier A carries the weight.** NCBI keyword searches for relB/relE return
~419,000 and ~304,000 records, and ~188,000 for both together. Almost all of
that is RelE/ParE superfamily background — chromosomes carry 5–12 copies each.
Sequence identity, not annotation, is what separates the lineage from the
background. This is also why the extended core anchor matters: the 527 bp seed
alone cannot discriminate, but `ISEsa1–ISEsa2–relBE` can.

**Nothing is deleted.** All 15 rows remain in `loci_selected.tsv` with
`selected` and `deselect_reason`. You can always see what the filter removed
and why. This is CLAUDE.md's "never silently drop" rule made mechanical.

**And the accounting is written out, not reconstructed.** Two intermediate
files make the whole cascade auditable:

- `blast/01_raw.tsv`, `blast/02_tier_a.tsv`, `blast/03_selected.tsv` — the
  search result at each stage, **all in the same `outfmt 6` columns**, so the
  narrowing can be measured with `wc -l` and inspected with `diff`:

```
  01_raw.tsv        15 HSPs   blastn, no filtering applied
  02_tier_a.tsv     14 HSPs   Tier A: pident>=95 qcov>=95; PARALOGUE removed
  03_selected.tsv    7 HSPs   Tier B: tier_a=FULL replicon_type=plasmid

  $ diff <(grep -v '^#' blast/01_raw.tsv | cut -f2,3) \
         <(grep -v '^#' blast/02_tier_a.tsv | cut -f2,3)
  < CP006642   67.681          <- the one hit sequence identity removed
```

  Each file carries its filter and its in/out counts as `#` comments, and
  `01_raw.tsv` additionally records the exact blastn command line — so a result
  can be checked against what BLAST actually returned rather than taken on
  trust.
- `filter_summary.tsv` — counts at every axis, each distinct combination of
  deselection reasons, and the parameters that produced them:

```
stage            category                   n   disposition
input            blastn HSPs                15
tier_a           FULL                       14  kept
tier_a           PARALOGUE                   1  dropped
replicon_type    chromosome                  7  dropped
replicon_type    plasmid                     8  kept
deselect_reason  replicon_type=chromosome    7  dropped
deselect_reason  tier_a=PARALOGUE            1  dropped
output           selected loci               7
params           tier_a=FULL replicon_type=plasmid min_flank=0 flank=20000
```

Note the filters are evaluated **together, not in sequence** — a hit can fail
on more than one — so the summary reports the population at each axis rather
than pretending to a strict funnel the code does not implement. Here 8 of 15
HSPs are on plasmids, but one of those is the paralogue, so 7 are selected.

Every threshold is a parameter — see §9.

---

## 4b. Design rule: no MGE feature is ever a membership criterion

**Lineage membership is decided by sequence identity. Context is recorded, never
used to admit or exclude a locus.** This holds for `czc`, for `Cytb/b6`, and
most importantly for ISEsa1/ISEsa2.

The tempting shortcut is "no vehicle found, therefore a paralogue". Measured on
the 16 curated records, that would be wrong in both directions that matter:

|  | has ISEsa1 | no ISEsa1 |
|---|---|---|
| **>=95% identity** (the lineage) | 12 | **2** |
| **<85% identity** (superfamily) | 0 | 33 |

Those 2 are `CP017220` (99.4%, class `NO CORE`, "IS1, adjacent deletions") and
`CP006642` (99.4%, carried by IS26 + Tn21 instead). Both are unambiguously the
lineage by sequence. Classifying them as paralogues would discard precisely the
evidence for **repeated independent acquisition** — the branch of the §6
decision tree that the project exists to resolve.

The vehicle is also far less reliable than a curated set suggests. Even among
nine loci *chosen* for showing the pattern:

```
COMPLETE   5/9    both ends intact
HALF-CORE  2/9    ISEsa2 absent entirely
DECAYED    1/9    ISEsa1 reduced to a 275 bp stub (17% coverage)
NO CORE    1/9    neither end
```

Only five of nine have an intact vehicle, and that is in a biased sample. In
`nt_prok` the fraction will be much lower — **what it actually is, is the
research question**, so it cannot also be the filter.

### Consequence for the extended core anchor

The core anchor (`ISEsa1–ISEsa2–relBE`) is built for **specificity**, and must
not be used for **discovery**:

| query | job | property |
|---|---|---|
| the 527 bp seed | find every copy | sensitive |
| the extended core | say which copies look ancestral | specific |

Searching with the core would miss both HALF-CORE loci, barely catch the
DECAYED one, and lose `CP017220` and `CP006642` entirely. It classifies what the
seed finds; it never decides what enters the pipeline.

**Enforced in code:** `select_loci.py` filters on `tier_a` (sequence) and
`replicon_type` (record) only. No filter anywhere consults IS or cargo content,
and a locus with zero IS hits still gets a row in every table.

**Demonstrated.** Running with `--tier_a all --replicon_type all` keeps every
hit. Three loci come back `NO CORE`, and only the sequence axis tells them
apart:

| locus | Tier A | identity | context |
|---|---|---|---|
| `CP006642__…__1` | FULL | 99.43% | NO CORE — carried by IS26 + Tn21 |
| `CP017220__…__1` | FULL | 99.43% | NO CORE — IS629, IS1R, ISSen4 present, no vehicle |
| `CP006642__…__2` | **PARALOGUE** | **67.68%** | NO CORE — a genuinely different gene |

Collapsing the axes would lump all three. Note also that `CP017220` is not
short of IS elements — it has three — they are simply not ISEsa1/ISEsa2. So
"has IS" is no more a valid lineage test than "has the vehicle" is.

> **Disabling a filter: use the literal `all`, not `''`.** Nextflow parses an
> empty CLI value as the boolean `true`, which then arrives as the string
> `"true"` and matches nothing — selecting zero loci with no error.

### Why Bakta does not replace the ISfinder layer

Annotation names IS elements by **family**, not by element. The same region of
CP028197:

| ISfinder (sequence homology) | GenBank / Bakta-style annotation |
|---|---|
| **ISEsa1** at −5494 | "ISNCY family transposase" |
| **ISEsa2** at −3898 | "IS110 family transposase" |

This is PIPELINE.md §0's warning made concrete: *"annotation keywords cannot be
trusted to find the vehicle — only sequence homology can."* Bakta gives the
complete gene track, boundaries, `cytb` and `czc`; ISfinder names the vehicle.
Both are needed, which is why the figure uses a two-tier colour key — ISfinder
element name where a feature overlaps a tiled IS call, gene family otherwise.

## 5. Flank availability, and why short flanks are kept

A hit near the end of a record cannot give a full ±20 kb window. The default
`min_flank = 0` **keeps** those loci and records how much was actually obtained.

```
circular, replicon big enough   → full flank both sides, window wraps the origin
circular, replicon too small    → whole replicon emitted once, whole_replicon=yes
linear                          → clamped at the ends; truncation is real
```

`flank_up`, `flank_dn`, `truncated_up`, `truncated_dn`, `wrapped` and
`whole_replicon` are all recorded, **post-reorientation** — they describe the
emitted sequence, not the replicon.

> **Unknown topology defaults to `linear`, deliberately.** BLAST output carries
> no topology field, so in a sweep it is usually unknown. Treating an unknown
> record as circular would *wrap* the window across the origin, splicing the
> end of the sequence onto its beginning and **fabricating a junction that does
> not exist**. Treating a circular record as linear merely clamps and records
> the lost flank. One error invents data; the other only loses context.
> `topology_evidence` marks which records were assumed rather than known.

> Truncation is itself a signal: an insertion hotspot breaks assemblies, and
> the *bla*NDM study lost 2,269 contigs exactly this way. Dropping short-flank
> loci discards that signal. It matters most for draft assemblies, where a hit
> lands on a contig that may itself be shorter than the window.

---

## 6. How flanking sequence is retrieved

Two mechanisms, **both implemented**. The coordinate handling is identical;
only the source differs, and the outputs are byte-identical either way.

**From source FASTAs** — `--fasta-dir`. Natural for curated sets where every
genome is already on disk.

**From the BLAST database** — `--blastdb`. The database *already contains the
full sequences*, so a hit's flanking DNA is already local:

```bash
blastdbcmd -db blastdb/public16 -entry CP018961.1 -range 237696-278222
# → 40,527 bp, byte-identical to the pipeline's window
```

> **This is not an optimisation, it is a requirement for pre-formatted
> databases.** `update_blastdb.pl nt_prok ...` ships no FASTAs at all, so
> `--fasta-dir` has nothing to read. `extract_loci.py` and `plasmidfinder.py`
> both fall back to `blastdbcmd` when a FASTA is unavailable.
>
> Verified: all 7 windows and all 17 replicon calls are identical whether
> sourced from FASTAs or from the database.

The database path is given as an **absolute path without file extensions**
(`--blastdb /scratch/$PROJECT/$USER/blastdb/nt_prok`) and is referenced in
place, never staged into the work directory — 127.7 GB cannot be copied per
task.

`efetch` by accession and range is needed in exactly one case: a hit found by
**remote** BLAST against a database you did not build. It needs network and is
rate-limited.

### The coordinate flow

```
blastn gives sstart/send on the subject
   → convert to replicon coordinates, normalise orientation
   → add ±flank, wrapping modulo replicon length if circular
   → clamp if linear, and record what could not be obtained
   → reverse-complement if the bait is on the minus strand
   → cut from the retrieved sequence
```

**Reorientation is not optional.** The bait sits on the minus strand in some
replicons and the plus strand in others. Without a common orientation every
downstream offset is mirrored for half the loci and nothing lines up in the
figure. All offsets assume the bait on `+`.

---

## 7. Outputs

| file | one row per | key columns |
|---|---|---|
| `results/bait_blast_raw.tsv` | BLAST HSP | the **verbatim** blastn output, with the exact command as a header comment — before any parsing or classification |
| `results/bait_hits.tsv` | BLAST HSP | `sacc, pident, qcov, strand, cov_*, tier_a, reason` |
| `results/filter_summary.tsv` | filter stage | `stage, category, n, disposition` — the cascade audited end to end |
| `results/record_metadata.tsv` | source record | `replicon_type, replicon_name, genus, slen, topology, is_complete, replicon_type_evidence` |
| `results/loci_selected.tsv` | hit | `locus_id, flank_up, flank_dn, truncated_*, wrapped, whole_replicon, selected, deselect_reason` |
| `results/loci.tsv` | extracted window | `window_len, bait_in_window, reoriented, source_interval` |
| `results/mge_hits.tsv` | IS call | `element, family, win_start, win_end, query_cov, cov_class, is_fragment, clipped_at_window_edge, clip_cause, clipped_missing_bp, offset_from_bait, nested_in, interrupts` |
| `results/replicons.tsv` | replicon hit | `accession, replicon, pident, qcovs, start, end` |
| `results/context.tsv` | locus | `vehicle_class, vehicle_ends, vehicle_interrupted, vehicle_edge_clipped, n_is_fragment, n_is_fragment_real, cargo_families, czc_offset, hok_gef_offset` |
| `results/blocks/*.blocks.tsv` | drawable feature | `feature_type, name, start, end, clipped_at_window_edge, colour_key, colour_source, annot_source` — the figure's input |
| `results/summary.tsv` | locus | one line per locus joining all of the above, incl. `lab_alias`, `plasmidfinder_replicons`, `has_ISEsa1`, `has_ISEsa2`, `nested_events` |
| `loci/plasmids/*.fa` | locus | the window itself; the FASTA header carries flank, truncation, wrap and bait position |

Read any of these with `scripts/tsvcols.py <file> <column> ...`, which selects
columns **by name**. `cut -f<numbers>` is not safe on these tables: inserting a
column silently shifts every later number, so a view that looks right returns
the wrong data. The same failure mode produced three real bugs in the scripts
themselves (see the `clipped_at_window_edge` commit), all now name-indexed.

`results/run_report.txt` is the place to start — one readable account of what
went in, what each filter removed, and what came out at every layer:
discovery and the narrowing, window statistics, the IS layer broken down by
coverage class, cargo families per locus, vehicle classes, replicon typing, and
the feature track with any locus that has no gene track called out by name.

`run_report.tsv` carries the same numbers for scripting, and `summary.tsv` is
the per-locus table.

`RUN_REPORT` reads only the **published** tables, so the report regenerates
from `results/` alone without re-running the pipeline:

```bash
python3 scripts/run_report.py --results results
```

Every emitted table has a header. (`results/stage1_is_hits.tsv`, from the
earlier hand-run work, is headerless raw BLAST in a non-standard column order
and is retained only as a historical artefact.)

---

## 8. Profiles and portability

| profile | tools from | for |
|---|---|---|
| `local` | **native arm64** binaries, **no containers** | the fast iterate-and-look loop |
| `firefly` | **Apptainer** + biocontainers | annotation, big-database layers |
| `gadi` | **Apptainer** + biocontainers, PBS Pro | the Stage 9 sweep only |

The local profile is deliberately asymmetric: biocontainers are linux-amd64 and
this machine's Docker daemon is off, so local runs use the native toolchain.

**No script contains a machine-specific path.** `scripts/toolpath.py` resolves
each executable in the order `--blast-bin` → `$BLAST_BIN` → `PATH`, and raises
an error naming all three if it fails. Local runs pass an explicit directory
because BLAST+ here is a native arm64 build deliberately kept off `PATH`;
cluster runs leave it empty and the binary comes from a container or module.

> `firefly.config` and `gadi.config` are **not yet verified against the real
> machines** — executor name, queue names and `-l storage` flags are
> best-guess. `env/firefly_probe.sh` collects what is needed to firm them up.

**Apptainer is the real reason to move to Linux**, more than CPU or internet.
`amrfinderplus` has *no macOS build at all* (bioconda or conda-forge, any
architecture), so Bakta is impossible locally; `keras` has no osx-arm64, so
geNomad is Rosetta-only. Both are fine in a linux-amd64 container.

**Internet is needed only at the edges** — one-time database pulls, the
metadata top-up (cacheable, `--offline_meta`), and the bulk sweep. Stages
between are fully offline, which is what makes Gadi's network-less compute
queues workable.

### Running on Gadi

Gadi provides everything needed, verified: `nextflow/24.04.5`,
`java/jdk-17.0.2`, `blast/2.11.0`, `singularity`, Python 3.13.

**Step 1 — one-time setup, on a LOGIN NODE** (it needs outbound network, which
standard compute queues do not have):

```bash
cd /g/data/rg42/TAflank
bash env/gadi_setup.sh
```

The repository deliberately does not carry the genome records
(`data/public_16/{fasta,gbk}`, 252 MB, recoverable by accession) or the two
reference databases, which are build artefacts. Without this step a run fails
at `BUILD_DB` with an empty input directory. The script is idempotent.

**Step 2 — submit:**

```bash
qsub env/gadi_run.sh                                    # plasmids, defaults
qsub -v REPLICON=chromosome env/gadi_run.sh             # chromosomes
qsub -v TIER=all,REPLICON=all env/gadi_run.sh           # keep everything
qsub -v DB=/scratch/rg42/blastdb/nt_prok env/gadi_run.sh  # the sweep
qsub -v BAKTA=true,BAKTA_DB=/g/data/rg42/bakta_db/db-light env/gadi_run.sh
```

Monitor with `qstat -u $USER`; output lands in `/scratch/rg42/$USER/taflank_run`.

### Two Gadi profiles, because scale changes the right shape

| profile | executor | for |
|---|---|---|
| `gadi_single` | `local` **inside one PBS job** | curated sets, tens of loci |
| `gadi` | `pbspro`, tasks fan out | the `nt_prok` sweep |

`gadi_single` is the default in `env/gadi_run.sh`. Everything runs inside the
job's own allocation, so there is no head process to keep alive on a login
node, no nested submission, and no head/task walltime coordination. It sizes
Nextflow's local executor from `PBS_NCPUS`/`PBS_VMEM` and puts the work
directory on `PBS_JOBFS`.

`gadi` is needed only when tasks must spread across nodes. Its head process
then has to sit on a login node (in tmux) or in `copyq`, because standard
compute queues have no outbound network — the same constraint that makes
`--meta_source blast` the default on both profiles.

### What changes on a cluster

Three things, each a parameter rather than a separate code path:

| local | cluster |
|---|---|
| builds its own BLAST db from FASTAs | `--blastdb <abs path>`; `BUILD_DB` is skipped |
| reads sequence from `--fasta-dir` | reads it from the database via `blastdbcmd` |
| `esummary` for record metadata | `--meta_source blast`, parsed from `stitle`, **no network** |

> **`-parse_seqids` cuts both ways.** It is load-bearing for genome databases —
> it is what lets `blastdbcmd` pull flanking sequence back out by accession and
> range, which is how `nt_prok` is used without shipping any FASTAs. But
> `makeblastdb` **refuses** it on the ISfinder mirror, which has 400 duplicated
> ids. That database is only ever a BLAST subject, never retrieved from, so it
> is built without.

### Databases: what to download and where to put it

Sizes are from NCBI's own metadata, not estimates:

| database | on disk | sequences | content |
|---|---|---|---|
| `nt_prok` | **127.7 GB** | 8.96 M | prokaryotic subset of nt, 503 Gbp |
| `nt_viruses` | **68.7 GB** | — | viral subset |
| `taxdb` | small | — | enables `staxids`/`sscinames` |

`nt_prok taxdb` is the recommended pair. `nt_viruses` costs 68.7 GB for little
benefit here — prophage context is better established with geNomad on the
extracted windows than by BLASTing against viral sequence.

> **The Gadi quota trap.** `$HOME` on Gadi is capped at **10 GiB per user**, and
> `df` does not show it — `df` reports the shared filesystem, which looks like
> tens of terabytes free. A download there fails with
> `curl (23) Failed writing body`, which reads like a network fault but is the
> write being refused. Check with `quota -s`, not `df`.

Where it can go, measured on this account:

| location | free | verdict |
|---|---|---|
| `$HOME` | 10 GiB cap | **no** — this is the trap above |
| `/g/data/rg42` | 54 GiB | **no** — too small for 127.7 GB |
| `/g/data/ei36` | 64 GiB | **no** — too small |
| `/scratch/rg42` | ~9 TiB | **yes** |
| `/scratch/ei36` | ~10 TiB | **yes** |

> **`/scratch` is purged** — files not accessed for ~100 days are deleted.
> Either re-touch the directory periodically or accept re-fetching. `/g/data`
> is persistent but no allocation here is large enough, so `/scratch` is the
> only option without a quota increase.

**How long it takes — measured on Gadi, not estimated.** GCP serves the
database files uncompressed, so `--decompress` is a no-op for `--source gcp`
and the transfer is the full on-disk size.

| streams | throughput | note |
|---|---|---|
| 1 | 17.8 MB/s | TCP-window-bound over a ~200 ms link, not bandwidth-bound |
| 4 | **75.6 MB/s** | 12 GB of real 3 GB volumes in 159 s |

At the 4-stream rate `nt_prok` lands in **~28 minutes**. Parallelism is what
matters here: a single stream is limited by latency, so `--num_threads 8` is
doing the real work. The 10 h walltime in the job script is generous.

**Verified result of an actual download** (`/scratch/rg42/blastdb`):

```
nt_prok     120 GB   8,962,499 sequences, 502.7 Gbp   ← blastdbcmd -info OK
nt_viruses   65 GB   (not needed -- see above)
taxdb       281 MB   all three files, so --taxonomy works
```

`env/gadi_fetch_blastdb.sh` runs the download as a **`copyq`** job rather than
on a login node: 127.7 GB is hours of transfer, and `copyq` is one of the few
Gadi queues with outbound network — standard compute queues have none, which is
also why `--meta_source blast` exists.

### Bait search at sweep scale

Two deviations are worth making deliberately when the database is `nt_prok`
rather than a curated set:

**Use `--blast_task megablast` for the bait search.** The `blastn` default
(word size 11) exists because megablast misses IS copies at 85-92% identity --
but that applies to the **IS search**, not the bait search, where only `FULL`
hits at >=95% identity are carried forward. Against 502.7 Gbp, word size 11 is
far more expensive for no gain in what is kept.

> Honest cost: megablast will not report the `<85%` superfamily background at
> all, so the `PARALOGUE` class disappears rather than being logged. At this
> scale that is arguably desirable -- logging millions of paralogues helps
> nobody -- but it **is** a deviation from "rejections are logged, never
> silent" and belongs in the methods, not slipped in silently.

**Scope before committing to the full chain.** Every `FULL` hit becomes a
window extraction, an IS search and a PlasmidFinder run. Run `BAIT_SEARCH`
alone first, look at the hit distribution, set the Tier A/B thresholds against
what actually came back, then re-run with `-resume` so the search is not
repeated.

> **The 527 bp seed is not specific enough for this database.** `nt_prok` is
> 8.96 M mostly-chromosomal sequences and chromosomes carry 5-12 relE-family
> genes each. Tier A at 95/95 does a lot of work, but the proper discriminator
> is the extended core anchor (`ISEsa1-ISEsa2-relBE`), which is **not yet
> built** -- see section 11. A seed sweep now is reconnaissance, not the real
> sweep.

### An already-present alternative on Gadi

`/g/data/xe2/references/blastdb_nt_April_2026` is a working full-`nt` database
(**125.7 M sequences, 3.70 Tbp**, April 2026), world-readable and verified with
`blastdbcmd -db nt -info`. Zero download. Three caveats, the last probably
decisive:

- it is full `nt`, **7× larger** than `nt_prok` (3.70 Tbp vs 503 Gbp), and with
  `-task blastn` that is a very different runtime;
- it includes all eukaryotic sequence, which is noise here and adds to the
  paralogue burden Tier A already works hard against;
- `/g/data` is mounted on compute nodes only when declared with
  `-l storage=gdata/xe2`, which NCI normally restricts to project members.
  It is readable from a **login node**; that does not mean a `qsub` job can see
  it. Verify before relying on it.

### Why `./tools/nf` and not bare `nextflow`

Three machine-specific traps, all encountered for real:

1. JDKs here are keg-only, so `java` is not on `PATH`.
2. The shell profile exports `JAVA_HOME=/opt/homebrew/opt/openjdk@17` — the
   **arm64** Homebrew prefix — but Homebrew here is **x86_64** at `/usr/local`,
   so that path does not exist. An inherited-but-broken `JAVA_HOME` is worse
   than none, so the wrapper validates it rather than defaulting around it.
3. **On arm64 the JVM's architecture leaks to child processes.** Running under
   an x86_64 JDK makes every task x86_64, and `/usr/bin/python3` then dies with
   `xcrun: unable to load libxcrun … need 'x86_64'` because CommandLineTools
   here is arm64-only. The wrapper probes each candidate JDK's `os.arch` and
   refuses a non-native one unless there is no alternative.

Related: **do not `brew install` the bioinformatics stack on this machine.**
Homebrew is x86_64/Rosetta here, so brew's BLAST would be a *downgrade* from
the native arm64 BLAST 2.17.0 that produced every verified result.

---

## 9. Parameters

All settable on the command line (`--flank 15000`) or in `nextflow.config`.

| parameter | default | meaning |
|---|---|---|
| `bait` | `anchors/relBE_seed.fa` | the bait sequence — any gene |
| `subregions` | `relB:1-240,relE:240-527` | named bait parts for orphan calls; `''` to disable |
| `genome_dir` / `db_name` | `data/public_16/fasta` / `public16` | build a db from these FASTAs |
| `blastdb` | `''` | absolute path to a **pre-formatted** db (no extensions). Set it and `BUILD_DB` is skipped and the db is never staged |
| `taxonomy` | false | also request `staxids`/`sscinames`; needs `taxdb` beside the db |
| `max_target_seqs` | 100000 | **not** a top-N filter — it prunes during search, so it must exceed the hits you expect |
| `blast_cpus` | 1 | threads for the bait search |
| `meta_source` | `ncbi` | `ncbi` \| `cache` \| `blast`. **`blast` makes no network call at all** |
| `meta_cache` | `data/metadata_cache` | per-accession JSON cache |
| `full_pident` / `full_qcov` | 95 / 95 | Tier A FULL threshold |
| `paralogue_max_pident` | 85 | below this is superfamily background |
| `blast_task` | `blastn` | **not** megablast locally — but see the sweep note below |
| `tier_a` | `FULL` | Tier A classes to carry forward |
| `replicon_type` | `plasmid` | `''` keeps all |
| `min_flank` | 0 | 0 keeps short-flank loci and records them |
| `flank` | 20000 | window half-width |
| `is_min_length` | **0** | asserted to stay 0 |
| `is_max_target` | 10000 | above the 5,970-sequence ISfinder db |
| `pf_min_cov` / `pf_min_id` | 80 / 95 | Acman et al. 2022 |


---

## 10. Validation

The pipeline is gated against a result verified by hand before it existed.

**Regression gate — CP018961 (pEC422_1).** All 7 rows of the archived
`mge_hits.tsv` are reproduced with **identical `offset_from_bait`**, while
`win_start`/`win_end` shift by exactly +5,000 — the expected consequence of
moving from ±15 kb to ±20 kb. Specifically:

- the **64 bp ISEsa1 fragment** is recovered at offset −6739;
- **ISKpn21 is called nested in ISEsa1**, matching slide 13.

**Negative controls.** CP006636 and CP018962 yield zero relBE, as documented.

**Paralogue separation.** The one RelE/ParE superfamily hit in CP006642
(67.7% identity, 90.6% of relE, 0.8% of relB) is classified `PARALOGUE` and
reported rather than silently carried.

**Window placement.** All 7 windows are 40,527 bp with the bait at exactly
20001–20527. The three windows whose bait differs from the seed differ by 3, 1
and 1 bases — exactly as their 99.43%, 99.81% and 99.81% identities predict.

**Independent agreement with the deck.** Slide 13 states that pECL_A is a
"variant FII FIB plasmid, own replicon types in PlasmidFinder". The pipeline
independently typed LT991956 as `IncFIB(pECLA), IncFII(pECLA)` from sequence
alone.

**Cluster mode equivalence.** The pipeline was re-run in the exact
configuration a cluster uses — pre-formatted database, `--genome_dir ''`,
`--meta_source blast`, no network — and compared against the FASTA + `esummary`
run:

- `mge_hits.tsv` — **identical** (45 IS calls)
- replicon typing — **identical** (17 hits)
- `loci.tsv` — 23 of 25 columns identical; the two that differ are `organism`
  (two-word binomial from the title vs the fuller `esummary` string) and
  `topology` (absent from BLAST output, so assumed — see §5)

**Cross-platform reproducibility, measured.** Job `182451448` on Gadi
reproduced the laptop run **exactly** — `mge_hits` 45/45, `context` 7/7,
`replicons` 17/17, with identical element names, `offset_from_bait` values,
vehicle classes and replicon calls. That holds despite the two runs differing
in nearly everything underneath:

| | laptop | Gadi |
|---|---|---|
| architecture | arm64 macOS | x86_64 Linux |
| BLAST | 2.17.0 | 2.11.0 |
| Python | 3.9.6 | 3.13.13 |
| metadata source | `esummary` (network) | `stitle` (offline) |
| executor | local, no containers | PBS job, `gadi_single` |

Cost: **0.38 SU**, 1 min 26 s walltime, 1.13 GB of 32 GB requested.

**Unit tests.** `/usr/bin/python3 -I tests/test_relbe_window.py` — 25
assertions covering circular wrap-around, the whole-replicon case, truncation
coordinate frames, and an end-to-end check against the real CP018961 record.

### Something the pipeline corrected

PlasmidFinder reassigned two plasmids that `data/public_16/accessions.tsv`
records as **IncFIB**:

| plasmid | annotated | typed |
|---|---|---|
| CP028197 pGMI14-002_1 | IncFIB | **IncHI2, IncHI2A**, IncC |
| CP028537 pQnrB4_020042 | IncFIB | **IncHI2, IncHI2A**, IncN |

This matters for the slide-7 hypothesis (0/37 IncHI2 plasmids with relBE from
*E. coli* vs 51/108 without), because it changes which plasmids sit in the
IncHI2 denominator. CP018961 additionally types as **pNDM-MAR**, a direct link
to the *bla*NDM literature the method is modelled on.

---

## 11. What is not built yet

| step | status | note |
|---|---|---|
| Direct repeats / TSDs (F6) | **not built** | would restore the `DR_seq` column the archived table had. Slide 13 already holds the ground truth: a 16 bp `AGGAGAAGAAGGCCAT` family, plus `TTAA` (ISKpn38) and `GAG` (ISKpn28) |
| Extended core anchor (Stage 2) | **not built** | `anchors/relBE_core.fa` does not exist. This is the gate before any expansion, and the reason a sweep cannot be keyword-led |
| Uniform annotation (Bakta) | **wired, off by default** | `--run_bakta`; needs a container and `--bakta_db`. The merge it feeds is already built and validated against a GenBank stand-in |
| Gene families (clustering) | **not built** | blastp + connected components is Roary's `-i 0.9 -s` and needs no new tool |
| Context classification (Tier C) | **BUILT** | `CLASSIFY_CONTEXT`; reproduces the hand classification on all four overlapping loci |
| Structural-variant splitting (Stage 6) | **not built** | Acman's method. Do **not** vendor the reference implementation: its released `track_variants.R` sets `step <- 100`, which cannot resolve the 64 bp fragment the decisive test depends on |
| Colour-block figure (Stage 7) | **not built** | needs Stage 3 for the gene track |
| ISEScan, IntegronFinder, mob_suite, geNomad | **not built** | all reachable via Apptainer on Firefly |
| NCBI sweep (Stage 9) | **not built** | gated on Stage 2; ~26,539 complete Enterobacterales assemblies |

### Current scope

Plasmids only, from the 16-record set — **7 loci**, which is 4 from the `nine`
set plus 3 from `extra`. The extras arrived free with the bait search, and two
of them (CP023893, LT991956) are the **pECL_A Class II** architecture that is
absent from the nine, so the working set spans all three architecture classes
rather than two.
