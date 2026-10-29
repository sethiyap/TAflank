# What every TAflank output means

One page per table. For each one: what it is, one row of real output, and what
every column holds. Everything here was read off an actual run of the 16
records, not written from the schema.

Read any table with column names rather than numbers:

```bash
python3 scripts/tsvcols.py <file>                      # list its columns
python3 scripts/tsvcols.py <file> <column> <column>    # show those columns
```

`cut -f<numbers>` is unsafe on these tables — adding a column silently shifts
every later number, so a view that looks right returns the wrong data.

---

## Start here: `run_report.txt`

The whole run on one page: how many hits went in, what each filter removed, and
what came out of every layer. If a number anywhere else looks wrong, this is
where you find out which stage it went wrong at. `run_report.tsv` is the same
content, machine-readable.

```bash
cat results/run_report.txt
```

---

## The deliverables

### `summary.tsv` — one row per relBE copy, everything joined

The table to hand someone who wants the answer without the workings.

| column | meaning |
|---|---|
| `locus_id` | `{accession}__{contig}__{copy}__{strand}` — the key every other table joins on. One row per copy, never merged; paralogous loci must stay separate. |
| `accession` / `lab_alias` | NCBI accession, and the name used in the lab or on the deck (`pEC422_1`) |
| `organism` | source organism as NCBI states it |
| `replicon_type` | `plasmid` or `chromosome` |
| `replicon_name` | the replicon's own name from the record (`pGMI14-002_1`) |
| `plasmidfinder_replicons` | Inc type(s) called by PlasmidFinder, e.g. `IncHI2`. This is **measured**, and it has corrected the deck: CP028197 and CP028537 are IncHI2, not IncFIB. |
| `replicon_len` / `topology` | full length of the parent replicon, and `circular`/`linear` |
| `window_len`, `flank_up`, `flank_dn` | how much sequence was actually obtained — not what was requested |
| `truncated_up` / `truncated_dn` | `yes` if the replicon ran out before the full flank. Recorded, never a reason to drop a locus. |
| `n_is` | IS elements found in the window |
| `n_is_complete` | of those, ≥90% of the reference element's length |
| `n_is_fragment` | <50% coverage |
| `n_is_fragment_clipped` | fragments that are only short because the window edge cut them — **not** evidence of decay |
| `n_is_fragment_real` | `n_is_fragment` minus the clipped ones. **This is the number the Stage 5 test counts.** |
| `is_elements` | comma-separated names, in position order |
| `has_ISEsa1` / `has_ISEsa2` | the two candidate vehicle ends |
| `nested_events` | `X>Y` means X sits inside Y — a later insertion that broke an earlier one |

### `mge_hits.tsv` — one row per IS element call

The MGE layer, every locus pooled. Per-locus copies are in `is/*.is.tsv` with
identical columns.

| column | meaning |
|---|---|
| `element` / `family` / `group` | ISfinder name, e.g. `ISEsa1` / `ISNCY` / `IS1202`. The **element** name is the vehicle; the family alone cannot identify it. |
| `win_start` / `win_end` | position in the window (1-based) |
| `length` / `strand` | aligned length, orientation |
| `pident` | % identity to the ISfinder reference |
| `query_cov` | % of the reference element covered |
| `cov_class` | `complete` ≥90%, `partial` 50–90%, `fragment` <50% |
| `is_fragment` | `yes` when `query_cov < 50`. A plain coverage statement — read it with the next three columns. |
| `clipped_at_window_edge` | `no` / `left` / `right` / `both` |
| `clip_cause` | **`-`** a real remnant in open DNA (evidence); **`window`** our ±20 kb cut sliced it (not evidence); **`replicon`** the replicon genuinely ends there (real truncation) |
| `clipped_missing_bp` | how much of the element lies beyond the edge |
| `offset_from_bait` | bp from the left edge of the bait. **Negative = upstream, positive = downstream.** This is the number the decisive test plots. |
| `sub_start` / `sub_end` / `slen` | which part of the reference element matched, and its full length |
| `nested_in` / `interrupts` | the element this one sits inside, if any |

**No length cutoff is applied anywhere in this layer.** The whole
ancestral-vehicle argument rests on a 64 bp ISEsa1 remnant, so a minimum length
would destroy the evidence. `is_min_length` is asserted to stay 0.

### `context.tsv` — one row per locus, the vehicle verdict

| column | meaning |
|---|---|
| `vehicle_class` | **COMPLETE** both ends present and intact · **DECAYED** both present, one broken — *the informative class* · **HALF-CORE** one end left · **NO CORE** neither; relBE arrived some other way |
| `vehicle_ends` | which of the ends were found |
| `vehicle_interrupted` | an end broken by a later insertion |
| `vehicle_edge_clipped` | an end whose true coverage is **unknown** because the window cut it. Such an end does not drive a DECAYED call; it means widen the window. |
| `n_is`, `n_is_fragment`, `n_is_fragment_real` | as in `summary.tsv` |
| `cargo_families` | which cargo families are present |
| `has_czc` / `has_cytb` / `has_hok_gef` | the genes the core is defined by, plus the second TA system |
| `czc_offset` / `hok_gef_offset` | their offsets. `czc` sits at **−1184** in every carrier — exactly, not approximately. |
| `metal_module` | `yes` when both `czc` and `cytb` are present |

Cargo presence is **recorded, never a membership test**. Requiring `czc` would
discard CP028197; `Cytb/b6` is absent from pEC422_1 altogether.

### `blocks/*.blocks.tsv` — the figure's input

One drawable feature per row: IS calls, annotated CDSs, and cargo hits on a
shared coordinate system.

| column | meaning |
|---|---|
| `feature_type` | `IS` / `CDS` / `cargo` |
| `name` / `family` | what to label it |
| `start` / `end` / `strand` | window coordinates |
| `offset_from_bait` | same convention as above — what to plot against |
| `cov_class`, `is_fragment`, `clipped_at_window_edge`, `nested_in` | carried through from the IS layer so the figure can mark a decayed element differently from one the window merely cut |
| `colour_key` | **what to colour by** |
| `colour_source` | `isfinder` or `annotation`. Two tiers: the ISfinder **element name** wins wherever a feature overlaps an IS call by ≥50%, the annotation name otherwise. Without this, gene-family clustering merges ISEsa1 with other ISNCY transposases and the figure misreports the vehicle. |
| `annot_source` | `bakta`, `genbank` or `none` — where the CDS calls came from |

### `replicons.tsv` — PlasmidFinder replicon typing

Run on the **whole parent replicon**, not the window: Inc replication genes are
rarely within 20 kb of relBE, so a window-only run returns almost nothing.

| column | meaning |
|---|---|
| `accession` | parent record |
| `replicon` | the Inc type, e.g. `IncHI2`, `IncFIB(K)` |
| `taxon_group` | which part of the PlasmidFinder database the hit came from, read off the hit's own header. The mirror holds 9 groups: `enterobacteriales` (159 sequences) plus the Gram-positive sets `Inc18`, `Rep1`–`Rep3`, `RepA_N`, `RepL`, `Rep_trans`, `NT_Rep`. |
| `pident` / `qcovs` | identity and coverage of the rep gene. Acman thresholds: coverage ≥80, identity ≥95. |
| `start` / `end` / `strand` | where on the replicon |

**Known gap.** The whole database is searched at once, so there is no
"unsupported taxon" branch: PlasmidFinder covers no Cyanobacteria,
Actinobacteria or *Helicobacter*, and a record from one of those simply returns
**zero rows** — currently indistinguishable from a plasmid that genuinely
carries no recognised rep gene. For the 16 Enterobacteriaceae records this does
not arise, but the TADB3 set spans those taxa, so the two cases need separating
before Stage 4.

---

## The workings — how the hits were narrowed

### `blast/01_raw.tsv` → `blast/02_tier_a.tsv`

Identical BLAST `outfmt 6` columns before and after the sequence filter, so you
can `diff` them and see precisely which matches were removed and why. `01_raw`
carries the exact `blastn` command as a header comment.

### `bait_hits.tsv` — every BLAST hit, classified

Standard BLAST columns (`pident`, `length`, `qstart`…`bitscore`, `qlen`, `slen`,
`stitle`) plus:

| column | meaning |
|---|---|
| `sacc` | subject accession. **Not `sseqid`** — with `-parse_seqids` that renders as `gb|CP018961.1|`. |
| `qcov` | % of the bait covered **by this one HSP** (there is no grouping by accession) |
| `strand`, `hit_start`, `hit_end` | orientation and position on the subject |
| `cov_relB` / `cov_relE` | % of each named bait part covered. Present only when `--subregions` is given; reporting only, it changes no decision. |
| `tier_a` | **FULL** ≥95% identity and ≥95% coverage · **PARALOGUE** <85% identity, the RelE/ParE superfamily background · **ORPHAN_\*** one named part survived · **PARTIAL** neither |
| `reason` | why, in words — every rejection is logged, never silent |

### `filter_summary.tsv` — the cascade audited

`stage`, `category`, `n`, `disposition`. Four columns that account for every
input hit through every filter, so no record disappears unexplained.

### `record_metadata.tsv` — what each hit sits on

NCBI `esummary` fields (`title`, `organism`, `taxid`, `slen`, `topology`,
`strain`, `biosample`…) plus TAflank's own calls:

| column | meaning |
|---|---|
| `replicon_type` | `plasmid` / `chromosome` / `unknown` |
| `replicon_type_evidence` | **how it was decided** — `esummary_genome` means NCBI stated it outright; `title_plasmid`, `title_complete_genome` etc. mean it was read off the description. Needed because CP021960, a 5.6 Mb circular chromosome, returns `genome: None`. |
| `is_complete` | complete record vs draft assembly |

### `loci_selected.tsv` — which hits became loci, and why not

Every Tier A hit with the Tier B decision attached: `selected` (`yes`/`no`) and
`deselect_reason`. Also `topology_evidence` — an unknown topology defaults to
**linear**, because wrapping a record whose topology you do not know fabricates
a junction that does not exist.

### `loci.tsv` — one row per extracted window

| column | meaning |
|---|---|
| `hit_start_orig` / `hit_end_orig` / `strand_orig` | where the bait was in the original record, before reorientation |
| `reoriented` | every window is flipped so relB is on `+` first. All offsets depend on this. |
| `window_len`, `flank_up`, `flank_dn` | what was actually obtained per side |
| `bait_in_window` | where the bait sits in the extracted window (`20001-20527` for a full ±20 kb) |
| `truncated_up` / `truncated_dn` | the requested flank was not achieved on that side |
| `wrapped` | the window crosses the origin of a circular replicon — handled with a modulo, so this is **not** truncation |
| `whole_replicon` | the window is the entire replicon, which is a better description than "truncated at both ends" |
| `source_interval` | true coordinates in the parent record |

### `cargo/*.cargo.tsv` — context genes

| column | meaning |
|---|---|
| `family` | one of 12 curated families (`cytb`, `czc`, `relB`, `relE`, `hok_gef`, `tra`, `qnr`, `intI`, `mer`, `ars`, `sil`, `ter`) |
| `source` | **the record the reference gene came from.** The most informative column in the file: for CP028197, `cytb`/`czc`/`relB`/`relE` are 100% identical to copies from four different records across three genera — which is what shows the block moved horizontally as one unit. |
| `pident` / `cov` | identity, and % of the reference gene covered |
| `clipped_at_window_edge`, `clip_cause`, `clipped_missing_bp` | as in the IS layer. A gene across the boundary reports low `cov` because we stopped looking, not because it is broken. |
| `offset_from_bait` | position relative to the bait |

This is deliberately **not** an annotation pass — 12 families answer the context
question in milliseconds. Full annotation is Bakta's job.

### `pf/*.pf.tsv`

Per-record PlasmidFinder hits, same columns as `replicons.tsv`.

### `bakta/*.gff3`

Full annotation of each window, when `BAKTA=true`. Optional, and complementary
to the IS layer rather than an alternative: Bakta names IS elements by **family
only** ("ISNCY family transposase"), so it cannot identify the vehicle —
ISfinder does that by homology. Both feed `blocks/`, where the element name wins
on ≥50% overlap.

Install it with `qsub env/gadi_bakta_setup.sh`; see RUNBOOK STEP 9b.

Without it, `annot_source` reads `genbank` and the CDS track comes from the
parent record's own annotation. That works for 15 of the 16 records; **LT991956
has 0 CDS features** and needs Bakta before it can contribute a gene track.

---

## Conventions that apply to every table

- **One record per copy.** Loci are named `{accession}__{contig}__{copy}__{strand}` and paralogous relE-family loci never collapse into one row.
- **Reorientation first.** Every window is flipped so relB is on `+` before anything else runs. All offsets assume it.
- **Negative offsets are upstream**, positive downstream, measured from the left edge of the bait.
- **Nothing is silently dropped.** Contig-truncated loci and fragmentary IS hits are kept and flagged. Truncation is itself a signal — the blaNDM study lost 2,269 contigs because the insertion hotspot breaks assemblies. Fragments are the evidence, not noise.
- **Every table has a header.**
