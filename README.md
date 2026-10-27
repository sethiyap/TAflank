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
                                                                          └─► collated tables
```

The method follows **Acman et al. 2022**, *Nat Commun* **13**:1131 —
[*Role of mobile genetic elements in the global dissemination of the carbapenem
resistance gene blaNDM*](https://doi.org/10.1038/s41467-022-28819-2).

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
    --blastdb /scratch/$PROJECT/$USER/blastdb/nt_prok \
    --genome_dir '' --meta_source blast --taxonomy
```

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
| **`results/summary.tsv`** | **locus — read this one first** |

## Requirements

- Nextflow ≥ 23.04 and a **native** JDK 17
- BLAST+ ≥ 2.17
- Python 3 (standard library only — no biopython, no pandas)
- R only for the figures, which are not built yet

Databases: an ISfinder mirror (`db/isfinder/`) and the CGE PlasmidFinder
database (`db/plasmidfinder/`), both with `PROVENANCE.txt` recording source and
fetch date.

## Documentation

| file | holds |
|---|---|
| `PIPELINE_DESCRIPTION.md` | **the implementation** — what runs, with what parameters, producing what files |
| `PIPELINE.md` | **the scientific rationale** — why each filter exists, the MGE layers, the decision tree |
| `CLAUDE.md` | project plan and conventions |

## Status

Working end to end on a curated plasmid set, gated against a result verified by
hand before the pipeline existed: all 7 archived IS calls for CP018961 are
reproduced with identical offsets, including the 64 bp ISEsa1 fragment and the
ISKpn21 nesting.

Not yet built: direct-repeat scanning, the extended core anchor, uniform
annotation and gene families, structural-variant splitting, and the
colour-block figure. See §11 of `PIPELINE_DESCRIPTION.md` — the gaps are listed
explicitly rather than implied by their absence.

## Citing

If you use TAflank, please cite Acman et al. 2022 for the method it follows.
The structural-variant splitting step is a reimplementation, not a fork, of
[`macman123/track_structural_variants`](https://github.com/macman123/track_structural_variants) —
see `PIPELINE_DESCRIPTION.md` for why, and for the parameter differences.
