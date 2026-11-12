# TAflank

**Bait-anchored discovery of mobile genetic elements in gene flanking regions.**

Give TAflank a gene — the *bait* — and it finds every copy of it across a set of
genomes, works out whether each copy sits on a plasmid or a chromosome, cuts out
the neighbourhood around it, and characterises the mobile genetic elements in
that neighbourhood.

It was built to answer a specific question about the *relBE* toxin–antitoxin
system, but **nothing in it hardcodes relBE**: point `--bait` at any gene and
the whole chain works.

```
bait.fa ──► local BLAST ──► Tier A ──► fetch metadata ──► Tier B ──► ±20 kb windows
 (any gene)   (vs a DB)    (sequence)    (batched)       (record)         │
                                                                          ├─► IS elements (ISfinder)
                                                                          ├─► replicon type (PlasmidFinder)
                                                                          ├─► gene annotation (Bakta / GenBank)
                                                                          ├─► pangraph: N-way alignment graph
                                                                          │   over all windows (sequence homology,
                                                                          │   independent of gene calls)
                                                                          └─► collated tables + figures
```

The method follows **Acman et al. 2022**, *Nat Commun* **13**:1131 —
[*Role of mobile genetic elements in the global dissemination of the carbapenem
resistance gene blaNDM*](https://doi.org/10.1038/s41467-022-28819-2) — for the
overall anchored, alignment-based approach, and uses **pangraph**
([Noll, Molari & Neher 2023](https://doi.org/10.1099/mgen.0.001034),
*Microbial Genomics* **9**:001034) as the sequence-homology engine itself: a
true N-way alignment graph rather than the pairwise `dc-megablast` Acman's
own pipeline uses.

---

## Why this shape

The hard part is not finding the bait. It is **separating genuine copies from
the homologous background**, and then keeping every fragment of evidence that
survives.

- **Filtering is tiered and explicit.** NCBI keyword searches for `relB`/`relE`
  return ~188,000 records, almost all RelE/ParE superfamily background —
  chromosomes carry 5–12 copies each. Sequence identity, not annotation, is
  what discriminates.
- **Nothing is silently dropped.** Every hit is emitted with a `selected` flag
  and a `deselect_reason`. Truncated loci are kept and the actual flank lengths
  recorded, because truncation is itself a signal: an insertion hotspot breaks
  assemblies, and the blaNDM study lost 2,269 contigs exactly that way.
- **No length cutoff in the IS layer.** A ~64 bp remnant can be the decisive
  evidence for an ancestral vehicle, mirroring the ~75 bp ISAba125 remnant that
  identified Tn125 for blaNDM. `is_min_length` defaults to 0 and is *asserted*
  to stay there.

## Quick start

```bash
./tools/nf run main.nf -profile local
```

Against a pre-formatted database on a cluster, with no network at all:

```bash
./tools/nf run main.nf -profile gadi \
    --bait anchors/relBE_seed.fa \
    --blastdb /scratch/$PROJECT/$USER/blastdb/nt_prok \
    --genome_dir '' --meta_source blast --taxonomy
```

`--bait` defaults to `anchors/relBE_seed.fa`, so it's optional here — shown
explicitly because it's the one flag that changes for a different gene;
everything else on the line is sweep-scale plumbing (database path, metadata
source, taxonomy).

Use `./tools/nf`, not bare `nextflow` — the wrapper resolves a **native** JDK,
which matters because on arm64 the JVM's architecture leaks to child processes.

## Outputs

| file | one row per |
|---|---|
| `results/bait_hits.tsv` | BLAST HSP, with its Tier A class and reason |
| `results/record_metadata.tsv` | source record, with replicon type and the evidence for it |
| `results/loci_selected.tsv` | hit, with flank availability and why it was or was not selected |
| `results/loci.tsv` | extracted window |
| `results/mge_hits.tsv` | IS call, including fragments and nesting |
| `results/replicons.tsv` | PlasmidFinder replicon hit |
| `results/context.tsv` | locus, with **two independent verdicts**: `vehicle_class` (VEHICLE-INTACT/-DECAYED/-HALF/-ABSENT, from ISEsa1/ISEsa2 coverage only) and `core_class` (CORE-INTACT/-DEGRADED/-PARTIAL/-ABSENT, from cytb/czc/relB/relE/hok_gef — never the vehicle). The two can and do disagree: CP010512 is VEHICLE-HALF with a fully intact core |
| `results/blocks/*.blocks.tsv` | drawable feature — the figure's input |
| `results/gene_families.tsv` | CDS, clustered by sequence (90% id / 80% cov) so the figure can colour genes by homology instead of by product name |
| `results/structural_groups.tsv` | locus, grouped with others sharing ≥90% of the vehicle + conserved-core genes — a fast triage for "do these plausibly share one history", not the decisive sequence-level test |
| `results/junctions.tsv` | IS call, with its target-site duplication where one boundary pair was clean enough to search, and whether that TSD is shared with another locus |
| `results/pangraph/graph.json` | pangraph's own N-way alignment graph (Noll, Molari & Neher 2023) over all loci — every input is an explicit ordered walk over blocks, so rearrangements are represented directly, not approximated |
| `results/pangraph_pairwise.tsv` | locus pair, with `shared_frac` (sequence identity-derived) — always produced, independent of `run_viz` |
| `results/pangraph_groups.tsv` | locus, grouped by `scripts/groups.py`'s union-find over `pangraph_pairwise.tsv` at `--min-similarity 0.85` |
| `results/pangraph_blocks.tsv`, `results/pangraph_context_groups.tsv` | the drawable pangraph blocks and the up/downstream context groupings behind the figures below |
| **`results/summary.tsv`** | **locus — read this one first** |

`figures/colour_blocks.png` (+ `.svg`) is the bait-aligned architecture panel —
one row per locus, IS and cargo coloured by role, an idealised-architecture
schematic computed from the data underneath. `figures/is_matrix.png` is the
companion IS presence/absence heatmap. Built by `scripts/colour_blocks.py`;
`--style tracks` switches from the default one-row-per-locus view to IS/CDS/
cargo on separate lines coloured by family.

`figures/pangraph_blocks.png` (+ `.svg`) is the sequence-homology counterpart —
one row per locus, coloured by pangraph block identity rather than gene name,
built by `scripts/pangraph_blocks.py` (`--row-scale`, default 0.5, packs rows
to roughly half height so more loci fit in one figure). `figures/
pangraph_annotated.png` overlays `colour_blocks.py`'s gene/IS arrows (Bakta's
calls among them, when `run_bakta=true`) on top of that same sequence-identity
background, one coordinate system, via `scripts/pangraph_annotated.py`.
`figures/pangraph_context_summary.png` (+ `.svg` and a two-sheet
`results/pangraph_context_summary.xlsx`) summarises loci by shared up/
downstream context group. This is now the **sole** sequence-homology method —
an earlier pairwise dc-megablast cross-check ("flanker") was retired once
pangraph's real alignment graph proved to subsume it; see `main.nf`'s comment
above `BUILD_PANGRAPH_GRAPH` for why.

`results/run_report.txt` accounts for the whole run on one page.
**[`OUTPUTS.md`](OUTPUTS.md) describes every table and every column.**

Read these with `scripts/tsvcols.py <file> [column ...]`, which selects by name
and lists the available columns when given none. `cut -f<numbers>` is unsafe
here: adding a column shifts every later number, so the view silently returns
different data.

## Requirements

- Nextflow ≥ 23.04 and a **native** JDK 17
- BLAST+ ≥ 2.17
- Python 3 (standard library only — no biopython, no pandas) plus `matplotlib`
  and `scipy` for the figures (`colour_blocks.py`, `structural_groups.py`,
  `pangraph_blocks.py`, `pangraph_annotated.py`), plus `openpyxl` for
  `pangraph_context_summary.py`'s workbook
- `pangraph` (native binary, no build needed — `env/gadi_pangraph_setup.sh`
  / `env/firefly_pangraph_setup.sh` fetch it) for the sequence-homology
  alignment graph
- Bakta (`db-light`), for uniform annotation where GenBank's own is absent or
  inconsistent — runs on Gadi via Apptainer (`env/gadi_bakta_setup.sh`);
  compute nodes have no network, so the container and database are pulled on
  `copyq` first and passed to the run as a local `.sif` path, never a
  `docker://` URI
- Same container-on-copyq pattern for the matplotlib/scipy stage on Gadi
  (`env/gadi_viz_setup.sh`) — see `env/gadi_run.sh`'s `VIZ=true` flag.
  Without it, the pipeline still produces every TSV table, just not the two
  `.png`/`.svg` figures from the sequence-homology stage.
- R is **not** required — the figures are Python/matplotlib, not R

Databases: an ISfinder mirror (`db/isfinder/`) and the CGE PlasmidFinder
database (`db/plasmidfinder/`), both with `PROVENANCE.txt` recording source and
fetch date.

## Documentation

| file | holds |
|---|---|
| `RUNBOOK.md` | **how to run it on Gadi**, one copy-paste step at a time, with full paths |
| `OUTPUTS.md` | **what every output table and column means** |
| `PIPELINE_DESCRIPTION.md` | **the implementation** — what runs, with what parameters, producing what files |
| `PIPELINE.md` | **the scientific rationale** — why each filter exists, the MGE layers, the decision tree |

`CLAUDE.md` (project plan and conventions) is kept locally and gitignored —
not published to GitHub.

## Status

Working end to end on **all 14 loci with a bait hit** (7 plasmid + 7
chromosome; `REPLICON=all`), gated against a result verified by hand before
the pipeline existed: all 7 archived IS calls for CP018961 are reproduced with
identical offsets, including the 64 bp ISEsa1 fragment and the ISKpn21
nesting. Adding the 7 chromosomal loci left the 7 plasmid rows in
`context.tsv` byte-identical — the gate held a second time.

**Built since that gate:** uniform Bakta annotation (`db-light`, all 14 loci);
gene-family clustering (`gene_families.py`); the architecture figure and IS
heatmap (`colour_blocks.py`, two styles); direct-repeat/TSD scanning
(`relbe_junctions.py` — proven working on real data, inconclusive for the
vehicle itself, see `PIPELINE.md` §3 Layer 8); same-history grouping by core
architecture (`structural_groups.py`); the independent `core_class` axis
alongside `vehicle_class`, so "vehicle intact but core missing" and the
reverse are both nameable, not just visually inferable; and pangraph
(`build_pangraph_graph.py`/`pangraph_blocks.py`/`pangraph_annotated.py`/
`pangraph_context_summary.py`) — a true N-way sequence-homology alignment
graph, giving a second, independent grouping and three figures alongside the
gene-annotation ones.

**[`REPORT.md`](REPORT.md)** walks through all of this on the actual
16-record pilot: the full pipeline flowchart, the Tier A/B/C filter numbers,
every annotation layer, the Jaccard and pangraph groupings side by side
(§7–7b), and the figures themselves.

**Not yet built:** the extended core anchor (`anchors/relBE_core.fa` — Stage
2, still the pivot: the current `core_class`/`structural_groups.py` work
*around* this gap, they don't close it), the TADB3 expansion (Stage 4), the
decisive sequence-level structural-variant test (all-vs-all `dc-megablast`,
threshold-stepped, `ggtree` — Stage 5; the figure's schematic panel and the
Jaccard grouping are useful proxies, not substitutes, and say so in their own
docs), and species/backbone clustering (Jaccard + PERMANOVA via `vegan` —
Stage 6). See `PIPELINE.md` §8 for the full per-stage table.

## Citing

See [`CITATION.cff`](CITATION.cff) for how to cite TAflank itself. It cites
pangraph (Noll, Molari & Neher 2023) and Acman et al. 2022 **separately**, in
its `references:` section — TAflank wraps and orchestrates both, but created
neither, and the file deliberately does not fold them into TAflank's own
authorship. The structural-variant splitting step is a reimplementation, not
a fork, of
[`macman123/track_structural_variants`](https://github.com/macman123/track_structural_variants) —
see `PIPELINE_DESCRIPTION.md` for why, and for the parameter differences.
