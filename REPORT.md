# relBE and its mobile genetic element context — pilot report

**Project:** TAflank · **Owner:** Pooja Sethiya · **Pilot date:** 2026-10-26 to 2026-11-03
**Source deck:** `relBE_schematic.pptx` (19 slides, 19 Aug 2026); hand-drawn predecessor:
`RelBE_from_Sept2018_lab_meeting.pptx`

---

## 1. What is relBE

`relBE` is a type II **toxin–antitoxin (TA) system** — a pair of adjacent, co-transcribed
genes:

- **relB** — the antitoxin. A small (79 aa) ribbon-helix-helix (RHH) DNA/protein-binding
  protein. On its own it is unstable and degrades quickly inside the cell.
- **relE** — the toxin. A slightly larger (95 aa) ribosome-dependent endoribonuclease.
  Left unchecked, it cleaves mRNA positioned in the ribosomal A site and halts
  translation.

The two genes overlap by one base pair and are transcribed as a single operon — the
classic co-transcribed TA arrangement. Under normal growth, RelB binds and neutralises
RelE. If the cell loses the DNA element carrying `relBE` (a plasmid cured, a chromosomal
region deleted), the short-lived RelB decays faster than the stable RelE, RelE is
released, and the cell's own translation machinery is shut down. This is a **post-
segregational killing / addiction module**: the TA pair does not do anything useful on
its own, but it makes whatever DNA molecule carries it much harder to lose, because
losing it is lethal.

That addictive property is the reason this project cares about `relBE` at all. It is
**not itself a mobile genetic element** — it has no transposase, no recombinase, nothing
that moves DNA. It is *cargo*: a hitchhiker that rides on mobile elements because those
elements benefit from being harder to lose. The question this pilot asks is what it
rides on, and whether the same vehicle carried it into every locus where it now sits, or
whether it got picked up independently, more than once.

### 1a. General schematic: vehicle, cargo, and the conserved core

Every locus in this pilot has the same basic shape, read left to right (upstream →
downstream). The terms below are used throughout the rest of this report:

```mermaid
flowchart LR
    vehicle["VEHICLE\n(e.g. ISEsa1-like,\noften decayed)"]:::vehicleStyle
    cargo["OTHER CARGO\n(variable)"]:::cargoStyle
    subgraph core["CONSERVED CORE (the anchor)"]
        direction LR
        cytb["cytb/b6"] --> czc["czc"] --> relB["relB"] --> relE["relE"]
    end
    downstream["usually EMPTY"]:::emptyStyle

    vehicle --> cargo --> core --> downstream

    classDef vehicleStyle fill:#ffd6e8,stroke:#c2255c,color:#5c0f2e
    classDef cargoStyle fill:#d0bfff,stroke:#6741d9,color:#2f1a66
    classDef coreStyle fill:#ffe066,stroke:#e67700,color:#4d2e00
    classDef emptyStyle fill:#eceff4,stroke:#4c566a,color:#2e3440,stroke-dasharray: 4 3
    class cytb,czc,relB,relE coreStyle
```

- **The core** (`cytb`/`b6` → `czc` → `relB` → `relE`) is the one part present and
  conserved across essentially every locus in this pilot — it is the anchor every
  coordinate in this report is measured from.
- **The vehicle** is the mobile element thought to have actually carried `relBE` in —
  most often an ISEsa1-like insertion sequence, frequently found degraded or with a
  second IS inserted inside it (§5–§7b asks whether this is *one* ancestral vehicle,
  reused repeatedly, or several independent ones that happen to look similar).
- **Other cargo** is everything else riding alongside — other IS copies, integrons,
  resistance genes — which varies locus to locus and is not itself the question this
  pilot asks.
- **Downstream is usually empty.** Variation in this dataset is one-sided: almost every
  insertion sits upstream of the core, which is itself a finding, not an assumption (§4
  checks this explicitly rather than just not looking).

## 2. The pipeline (TAflank)

TAflank is a bait-anchored discovery pipeline, built around a single idea: **anchor
every comparison on a gene of known identity, not on an assembly's own annotation**.
Point it at a bait sequence and a set of genomes, and it:

1. finds every high-confidence copy of the bait by BLAST, discarding sequence-divergent
   paralogues before they can contaminate anything downstream;
2. works out whether each copy sits on a plasmid or a chromosome;
3. cuts out a defined window of flanking sequence around each copy, reorienting every
   window the same way so positions are comparable across loci;
4. characterises the mobile genetic elements and the other genes (cargo) found in that
   window, using several independent tools and filters, none of which discard evidence
   silently.

The method follows Acman *et al.* 2022 (*Nat Commun* **13**:1131), who used exactly this
anchored, alignment-based approach to show that the carbapenem resistance gene *blaNDM*
spread on one ancestral vehicle. This pilot repeats that logic for `relBE`.

### Pipeline flow

Colour marks the role each stage plays: grey = input, orange = the Tier A/B/C filter
(§3), blue = window extraction (§4), the four annotation passes of §5 each get their
own colour, gold = the two classification axes (§5), teal = the TSD/Jaccard history
tests (§6–7), sky blue = pangraph's independent sequence-homology check (§7b), dark
green = the deliverables (§8). GitHub/VS Code render this statically — no
click-through — but still support pan/zoom and the colour coding below.

**Where it ran:** everything in this pilot ran locally (arm64 macOS) except
**Bakta** (marked `[Gadi]` below), which has no macOS build and runs on Gadi
via Apptainer — compute nodes there have no outbound network, so the
container and database are pulled on `copyq` first (`env/gadi_bakta_setup.sh`).
The matplotlib/scipy figure stage (`colour_blocks.py`, the pangraph figures)
has the same container-on-copyq option available on Gadi
(`env/gadi_viz_setup.sh`) for cluster-scale runs, but ran locally here.

```mermaid
flowchart TD
    bait["527bp relBE bait\n+ 16 GenBank records"]:::input

    bait --> blast["blastn bait vs records\n15 HSPs"]:::tier
    blast --> tierA{"Tier A\nidentity/coverage filter"}:::tier
    tierA -->|"FULL — 14 kept"| window
    tierA -->|"PARALOGUE — 1 dropped\n(CP006642 background hit, 67.7% id)"| dropped["excluded\nnever reaches window extraction"]:::dropped

    window["reorient relB to '+'\ncut ±20kb window\nclamped + recorded at replicon ends"]:::window
    window --> isfinder
    window --> cargo
    window --> bakta
    window --> pangraphBuild

    tierA -->|"FULL — 14 kept"| replClass{"replicon type?\n(from assembly metadata)"}:::tier
    replClass -->|"plasmid — 7 loci"| replicon
    replClass -->|"chromosome — 7 loci"| noPF["skipped\nno Inc-type replicon to call"]:::dropped

    replicon["whole parent replicon\n(not just the window)"]:::input --> plasmidfinder

    subgraph annotate["§5 — four independent annotation passes"]
        isfinder["ISfinder blastn\nNO length/identity cutoff\n→ complete / partial / fragment"]:::isfinder
        cargo["cargo-gene blastn\nrelB/relE/cytb/czc/hok_gef/tra/ars/sil/intI/qnr\n(curated from this project's own 16 records)"]:::cargo
        plasmidfinder["PlasmidFinder blastn\nInc typing, cov≥80% id≥95%"]:::plasmidfinder
        bakta["Bakta db-light\nCDS annotation\n[Gadi -- Apptainer]"]:::bakta
    end

    bakta --> families["sequence clustering\n90% id / 80% cov\n→ gene families (not product names)"]:::bakta

    isfinder --> vehicle["vehicle_class\nVEHICLE-INTACT/DECAYED/HALF/ABSENT"]:::classify
    cargo --> core["core_class\nCORE-INTACT/DEGRADED/PARTIAL/ABSENT"]:::classify

    isfinder --> tsd["relbe_junctions.py\nDR/TSD scan on 74 IS calls"]:::history
    vehicle --> jaccard["Jaccard on vehicle+core signature\nunion-find @ 90%"]:::history
    core --> jaccard
    families --> jaccard

    tsd --> figure
    jaccard --> figure
    plasmidfinder --> figure
    figure["figures/colour_blocks.png\nfigures/is_matrix.png\nresults tables"]:::output

    pangraphBuild["pangraph build\nN-way alignment graph\n(Noll, Molari & Neher 2023)"]:::pangraph
    pangraphBuild --> pangraphGroups["pangraph_groups.tsv\nunion-find @ shared_frac ≥0.85"]:::pangraph
    pangraphBuild --> pangraphCtx["pangraph_context_summary\n11 unique contexts"]:::pangraph
    pangraphGroups --> pangraphOut
    pangraphCtx --> pangraphOut
    pangraphOut["pangraph_blocks.png\npangraph_annotated.png\npangraph_context_summary.png/.xlsx"]:::output

    classDef input fill:#eceff4,stroke:#4c566a,color:#2e3440,stroke-width:1px
    classDef tier fill:#ffd8a8,stroke:#e8590c,color:#7a3700,stroke-width:1px
    classDef dropped fill:#ffc9c9,stroke:#c92a2a,color:#7a0000,stroke-width:1px,stroke-dasharray: 4 3
    classDef window fill:#a5d8ff,stroke:#1864ab,color:#0a3d62,stroke-width:1px
    classDef isfinder fill:#b2f2bb,stroke:#2b8a3e,color:#0b4a1e,stroke-width:1px
    classDef cargo fill:#d0bfff,stroke:#6741d9,color:#2f1a66,stroke-width:1px
    classDef plasmidfinder fill:#99e9f2,stroke:#0c8599,color:#063c45,stroke-width:1px
    classDef bakta fill:#ffd6e8,stroke:#c2255c,color:#5c0f2e,stroke-width:1px
    classDef classify fill:#ffe066,stroke:#e67700,color:#4d2e00,stroke-width:1px
    classDef history fill:#63e6e2,stroke:#0b7285,color:#063539,stroke-width:1px
    classDef pangraph fill:#d0ebff,stroke:#1971c2,color:#0b3d66,stroke-width:1px
    classDef output fill:#63e6be,stroke:#087f5b,color:#063d2e,stroke-width:2px
```

### Reproducing this pilot yourself

Run from the pipeline's own directory — every path in `nextflow.config` is
relative to it (`cd TAflank` after cloning, not from wherever your shell
happens to start):

```bash
cd TAflank
./tools/nf run main.nf -profile local \
    --bait anchors/relBE_seed.fa \
    --replicon_type all \
    --outdir test_relBE_run/results --figdir test_relBE_run/figures
```

`--replicon_type all` is what makes this the full 14-locus pilot (7 plasmid
+ 7 chromosome) rather than the config's own default of `plasmid`-only (7
loci). Sending the output to `test_relBE_run/` rather than the default
`results/`/`figures/` keeps it separate from anything already in the repo.
See the README's **Installation** section first if any of the tools
(`blastn`, `pangraph`, the Python figure dependencies) aren't set up yet.

**Skipping the rebuild, with a caveat.** The command above lets `BUILD_DB`
rebuild `blastdb/public16` from `data/public_16/fasta` on every run — fast
for 16 genomes, and one less thing to keep in sync. To search an
**already-built** database instead (the point of `--blastdb` for a real
sweep against something like NCBI's `nt_prok`), add it and clear
`genome_dir` so `BUILD_DB` is skipped entirely:

```bash
./tools/nf run main.nf -profile local \
    --bait anchors/relBE_seed.fa \
    --replicon_type all \
    --blastdb "$(pwd)/blastdb/public16" --genome_dir '' \
    --outdir test_relBE_run/results --figdir test_relBE_run/figures
```

`--blastdb` takes an absolute path with **no file extension** — this is a
prefix for `public16.n*`, not a directory. The caveat: setting `--blastdb`
means `BUILD_DB` genuinely never runs, so if `blastdb/public16.*` doesn't
already exist (e.g. it was cleaned up as a regenerable build artifact, which
it is — see `.gitignore`), `blastn` fails with *"No alias or index file
found for nucleotide database"* rather than building one for you. Rebuild it
by hand first if needed:

```bash
./scripts/build_blastdb.sh data/public_16/fasta public16 'public16 (manual rebuild)'
```

## 3. The pilot test: 16 input records

The pilot used **16 closed GenBank records** as the search database — nine loci and two
negative controls hand-picked from a September 2018 lab-meeting survey, plus five
"extra" records that turned up free alongside them. Every record is annotated below with
how it was actually used, not just what it was labelled.

| accession | organism | set | replicon | relBE bait hit? |
|---|---|---|---|---|
| CP028197.1 | *Salmonella enterica* Concord | nine | plasmid pGMI14 (IncFIB) | yes, 100.0% id |
| CP028537.1 | *Enterobacter hormaechei* | nine | plasmid pQnrB4_020042 (IncFIB) | yes, 100.0% id |
| MF344582.1 | *Citrobacter freundii* NY2385 | nine | plasmid pNY2385-HI2 | yes, 100.0% id |
| CP018961.1 | *Escherichia coli* Ecol_422 | nine | plasmid pEC422_1 | yes, 100.0% id |
| CP010512.1 | *Enterobacter cloacae* colR/S | nine | **chromosome** | yes, 100.0% id |
| CP026216.1 | *Citrobacter* sp. CFNIH10 | nine | **chromosome** | yes, 100.0% id |
| CP020841.1 | *Klebsiella pneumoniae* KPN1482 | nine | **chromosome** | yes, 100.0% id |
| CP012599.1 | *Salmonella enterica* Newport | nine | **chromosome** | yes, 100.0% id |
| CP017220.1 | *Escherichia coli* FAM21845 | nine | **chromosome** | yes, 99.4% id |
| CP021960.1 | *Klebsiella pneumoniae* AR_0139 | extra | **chromosome** | yes, 100.0% id |
| CP026975.1 | *Enterobacter cloacae* FDAARGOS_77 | extra | **chromosome** | yes, 99.8% id |
| CP023893.1 | *Klebsiella ornithinolytica* FDAARGOS_431 | extra | plasmid (unnamed1) | yes, 99.8% id |
| LT991956.1 | *Enterobacter hormaechei* steigerwaltii C309 | extra | plasmid | yes, 99.8% id |
| CP006642.1 | *Escherichia coli* PCN061 | extra | plasmid PCN061p6 | yes, 99.4% id |
| CP018962.1 | *Escherichia coli* Ecol_422 (chromosome) | **negative control** | chromosome | **no — zero hits** |
| CP006636.1 | *Escherichia coli* PCN061 (chromosome) | **negative control** | chromosome | **no — zero hits** |

**14 of the 16 records carry relBE; 2 do not, by design.** CP018962 and CP006636 are the
*chromosomes* of two strains whose *plasmids* (CP018961 and CP006642 respectively) are in
the table above — they were included specifically to confirm the search does not call a
hit where there is none. Both returned **zero** BLAST hits against the 527 bp bait,
clean negatives.

**Replicon split of the 14 positive loci: 7 plasmid, 7 chromosome.**

### Was there a false positive?

One. BLAST found a **second**, much weaker hit inside CP006642 — 67.7% identity over
263 bp (half the bait length), sitting ~400 bp from the genuine 99.4%-identical,
full-length copy. This is a diverged paralogue: a different, chromosomally-background
`relE`-family gene (the toxin-antitoxin superfamily this bait belongs to is large —
keyword searches for `relB`/`relE` alone return on the order of 188,000 GenBank records,
almost all unrelated background). The pipeline's **Tier A** sequence-identity filter
exists precisely to catch this: hits are kept as `FULL` only above a high identity/
coverage threshold and classed `PARALOGUE` and dropped otherwise. The CP006642
paralogue was correctly classed `PARALOGUE` and excluded; only the genuine copy
proceeded to window extraction.

```
input      blastn HSPs                 15
tier_a     FULL (kept)                 14
tier_a     PARALOGUE (dropped)          1   ← the CP006642 background hit
```

## 4. How the flanks were determined

Once a locus passed Tier A, a fixed **±20 kb window** was cut around it with these
rules, applied identically to every locus:

- **Reorient first.** Every window is flipped, if needed, so `relB` sits on the `+`
  strand. All downstream offsets are measured in this orientation, so a position
  reported as "−1184 bp from relB" means the same physical relationship in every locus,
  regardless of which strand the gene was on in the original record.
- **Coordinate zero is the relB start codon**, not the edge of the window.
- **Clamped, not padded, at replicon ends.** A window near the end of a plasmid or
  contig gets however much real sequence exists there — never padded with Ns, and the
  shortfall is recorded explicitly, because an assembly breaking near this locus is
  itself a signal (insertion hotspots are exactly where assemblies tend to fragment),
  not noise to be hidden.
- **One record per copy**, named `{accession}__{contig}__{copy}__{strand}`, so that if a
  genome ever carries more than one relBE copy they are never silently merged into one
  node.

Nothing about this stage filters anything out — every locus that passed Tier A gets a
window, whatever its contents turn out to be.

## 5. Annotating what's in the window: tools and filtering criteria

Four largely independent passes were run over each window, with results kept
deliberately separate rather than collapsed into one score:

1. **Insertion sequences (IS)** — `blastn` vs a local ISfinder mirror, no length or
   identity cutoff.
2. **Cargo genes** — `blastn` against a small curated reference set drawn from this
   project's own 16 records.
3. **Replicon typing** — `blastn` against the CGE PlasmidFinder database, run on the
   whole parent replicon (plasmid loci only).
4. **Uniform gene annotation** — Bakta (`db-light`), re-annotating CDS from scratch
   across every window.

(Gene-family clustering, below, is a follow-on analysis built from pass 4's output,
not a fifth independent pass.)

### Insertion sequences (IS) — `blastn` vs a local ISfinder mirror

Every window was searched against the ISfinder database with **no length or identity
cutoff**. This is a deliberate departure from how IS-finding tools are normally run: a
64 bp fragment of a decayed element can be the single most important piece of evidence
in the whole analysis (see §6), and a default minimum-length filter would discard it
silently. Instead:

- hits are tiled (overlapping HSPs merged, 30 bp tolerance for shared bases at
  junctions) and checked for **nesting** — one IS landing inside another, which is
  itself a later, separate insertion event, not noise;
- every hit is classed by how much of the reference element it covers:
  **complete** (≥90%), **partial** (50–90%), or **fragment** (<50%) — and **nothing is
  filtered by this class**, fragments are kept as evidence, not discarded as noise.

Across the 14 loci this found **74 IS calls**: 56 complete, 2 partial, 16 fragments,
spanning **12 ISfinder families** (ISNCY ×15, IS3 ×14, IS6 ×11, IS110 ×9, IS481 ×6, and
smaller families).

### Cargo genes — `blastn` against a small curated reference set

A second, separate search targets a short list of named genes of specific interest to
this project, not an unbiased scan: `relB`, `relE`, `cytb` (cytochrome b561), `czc`
(cobalt–zinc–cadmium efflux), `hok_gef` (a second small toxin), `tra` (conjugal
transfer), `ars` (arsenic resistance), `sil` (silver/copper efflux), `intI`
(integrase), `qnr` (fluoroquinolone resistance). The reference sequences for each were
extracted directly from this project's own 16 records (traceable to an accession and
coordinates), not pulled from an external curated database.

**How the set was actually built** (`scripts/build_cargo_set.py`): for each family, a
hand-written regex is matched against every CDS's `/product=` text across all 16
`.gbk` files (e.g. `czc` matches `cobalt-zinc-cadmium|czc[ABCD]|cation efflux|CzcD|
cusA|cusB`). Every matching CDS (≥150 bp) is sliced straight out of that record's own
sequence, reverse-complemented if on the `complement` strand, and deduplicated by its
first 200 bp so near-identical alleles collapse to one representative while genuinely
distinct sequences under the same pattern are kept. Each entry is named
`{family}|{accession}|{locus_tag}|{coords}`, so it stays traceable to its source.
This is **product-name-regex curation, not manual sequence review** — it inherits
whatever each original GenBank submitter called the gene, so a cargo gene that is
unannotated or named unexpectedly in the source record is invisible to this set.

The current `assets/cargo_genes.fa` holds **333 sequences** and is actually wider than
the ten families named above: it also carries `mer` (mercury resistance, 30 entries)
and `ter` (tellurite resistance, 62 entries), picked up by the same build but not yet
mentioned in this report's prose. Per-family counts: `tra` 128, `ter` 62, `ars` 42,
`mer` 30, `relE` 18, `czc` 11, `hok_gef` 11, `sil` 9, `qnr` 7, `relB` 7, `intI` 5,
`cytb` 3.

### Replicon typing — `blastn` against the CGE PlasmidFinder database

Run on the **whole parent replicon**, not just the ±20 kb window (Inc-type replication
genes are rarely within 20 kb of relBE), at Acman's published thresholds (coverage
≥80%, identity ≥95%). Of the 7 plasmid loci, all 7 typed, returning 17 replicon hits
across Inc groups IncHI2A (×3), IncHI2 (×3), IncFIB(K) (×2), IncN (×2), IncX1 (×1) and
IncFIB(pNDM-Mar) (×1) — several plasmids carry more than one replicon marker.

### Uniform gene annotation — Bakta (`db-light`)

GenBank's own CDS annotation across these 16 records ranged from richly annotated to
**entirely absent** — one record (LT991956.1) arrived with zero CDS features. Rather
than work with whatever each submitter happened to provide, every window was
re-annotated from scratch with **Bakta**, using its lighter reference database
(sufficient for CDS coordinates and most product names; the full database would add
~80 GB for a modest further gain in naming, not in the coordinates anything downstream
actually depends on). Bakta runs via Apptainer on the Gadi HPC cluster — its
AMRFinderPlus dependency has no native macOS build, so it cannot run on the project
laptop at all.

Result: **715 CDS called across the 14 loci**; 264 (37%) returned only "hypothetical
protein" as a product name — expected for AT-rich, horizontally-mobile DNA, which is
exactly the kind of sequence under-represented in reference databases.

### Gene families — sequence clustering, not product names

Because over a third of CDS have no usable name, genes were additionally clustered by
**sequence** (90% identity, 80% mutual coverage — the Roary `-i 0.9 -s` convention),
giving every CDS a family ID independent of whether it has a name at all. This is what
lets the figure (§8) colour "hypothetical protein" genes consistently when the same
unnamed gene recurs across loci, rather than painting every unnamed gene the same
flat grey. 715 CDS resolved into **470 families**, 116 of them shared by more than one
locus; the 264 previously-identical-looking "hypothetical protein" CDS split into 187
distinct families once sequence, not the product string, was used to group them.

### Two independent classification axes, deliberately never merged

Each locus gets **two separate verdicts**, because collapsing them into one score would
hide real biology:

- **`vehicle_class`** — is the IS "vehicle" (the pair of insertion sequences bracketing
  the conserved region, named `ISEsa1`/`ISEsa2`) intact? `VEHICLE-INTACT` (both present,
  ≥90% coverage each) / `VEHICLE-DECAYED` (both present, one degraded) /
  `VEHICLE-HALF` (only one present) / `VEHICLE-ABSENT` (neither).
- **`core_class`** — is the conserved cargo cassette (`cytb`, `czc`, `relB`, `relE`,
  `hok_gef`) intact? `CORE-INTACT` / `CORE-DEGRADED` / `CORE-PARTIAL` / `CORE-ABSENT`,
  scored the same way but **never looking at the IS elements at all**.

These two axes genuinely disagree in this dataset: **CP010512** is `VEHICLE-HALF` (it
kept only one of the two flanking IS elements) yet fully `CORE-INTACT` — the conserved
cargo cassette survived even though half its vehicle did not. Across the 14 loci:
6 `VEHICLE-INTACT`, 1 `VEHICLE-DECAYED`, 5 `VEHICLE-HALF`, 2 `VEHICLE-ABSENT`; 8
`CORE-INTACT`, 6 `CORE-PARTIAL`, 0 fully absent on either axis.

## 6. Insertion sequences, cargo, and direct repeats — in brief

- **Insertion sequences (IS)** are the simplest autonomous mobile genetic elements: in
  general, a single transposase gene flanked by short terminal inverted repeats, with no
  other cargo of their own. ISfinder groups them into families (IS3, IS5, IS6, IS110,
  IS481, "ISNCY" = not yet classified into a family, and others) by the transposase's
  sequence and mechanism. Most families cut-and-paste or copy-and-paste themselves into
  a new site by a staggered double-strand break; a few families (notably **IS110** and
  **IS481**, both present in this dataset) instead use a site-specific recombinase
  mechanism that does not leave the classic signature described next.

- **Cargo** is whatever gets carried *by* a mobile element or *alongside* it in the same
  neighbourhood, without itself providing any mobilisation function. In this project,
  `relBE` is cargo; so are the metal-resistance genes (`cytb`, `czc`) found immediately
  upstream in most loci, and, more sporadically, antibiotic resistance (`qnr`) and
  conjugal transfer machinery (`tra`). Cargo is the payload; IS elements, transposons and
  plasmid backbones are the vehicle.

- **Direct repeats / target-site duplications (TSDs)** are the forensic signature a
  classic (staggered-cut) transposase leaves behind. When it inserts, it cuts the target
  DNA with a short staggered overhang; repair of that overhang duplicates a few base
  pairs of **host** sequence, so an identical short sequence ends up flanking *both*
  ends of the newly-inserted element. Because each independent insertion event cuts at
  its own, essentially random, position, **independent insertions of even the identical
  IS element almost always leave different TSD sequences** — while every descendant of
  one single ancestral insertion event inherits the *same* TSD, because it is simply
  copied along with the rest of the DNA. Comparing TSDs across loci is therefore a
  direct test of this project's central question: a TSD shared, letter for letter, by
  several loci is strong evidence of one shared historical event; TSDs that differ from
  locus to locus point to repeated, independent acquisition.

  This pilot ran that scan (`relbe_junctions.py`) on all 74 IS calls. It is
  proven to work — it recovered a **100%-identical, 3 bp TSD (`ACC`) at the identical
  window coordinates** for an `ISKpn8` copy shared between CP012599 and CP028197,
  strong evidence those two loci inherited that particular insertion from one shared
  event rather than acquiring it twice. For the vehicle elements themselves
  (`ISEsa1`/`ISEsa2`), the scan could not reach a conclusion: wherever both are present
  they sit directly adjacent to each other, with no host-DNA gap between them to
  search — a structural finding in its own right (the vehicle behaves as one tight
  compound unit, not two independently-landed pieces), but one that leaves the direct
  ancestral-vehicle question still open pending a boundary-precision tool (ISEScan) this
  pilot has not yet run.

## 7. A fast proxy for "do these loci share one history"

Short of the full sequence-level structural-variant test (planned, not yet built — see
§9), a lighter check was run: **Jaccard similarity restricted to the vehicle and the
five conserved cargo genes** (not the full ±20 kb window, which is dominated by
locus-specific insertions and would hide the signal). At a 90% similarity threshold,
union-find clustering on this restricted signature split the 14 loci into:

| group | loci | shared architecture |
|---|---|---|
| 6 loci | CP012599, CP018961, CP020841, CP021960, CP028197, CP028537 | full vehicle + full core |
| 3 loci | CP023893, CP026975, LT991956 | ISEsa1 only, `hok_gef`, `relB`, `relE` |
| 2 loci | CP006642, CP017220 | bare `relB`/`relE` pair only |
| 3 singletons | CP010512, CP026216, MF344582 | each a genuinely distinct partial architecture |

Every multi-locus group above is **pairwise 100% identical** internally on this
restricted signature, not merely above the 90% cutoff — and notably, every one of the
three multi-locus groups **mixes plasmid and chromosomal loci**, consistent with
`relBE`'s vehicle moving between genomic contexts rather than architecture tracking
replicon type.

## 7b. A second, independent check: pangraph sequence homology

§7's Jaccard grouping asks "do these loci carry the same **named** elements" — it never
looks at the DNA itself, so two loci could share an element catalogue by convergence
rather than by inheritance. `pangraph` (Noll, Molari & Neher 2023, *Microbial Genomics*
**9**:001034) builds a true N-way alignment graph over all 14 windows and answers a
different, stricter question: how much of the **actual sequence**, not just which named
features, do two loci share.

At `shared_frac ≥ 0.85` over the **whole** ±20 kb window (not restricted to the
vehicle + core the way §7's signature is), pangraph finds far fewer multi-member
groups — these same groups are also the **11 unique up-/downstream contexts** drawn in
`pangraph_context_summary.png` (§8) and tabulated in
`results/pangraph_context_summary.xlsx`:

| context | n | representative (medoid) | members | weakest internal pair |
|---|---|---|---|---|
| 1 | 3 | CP028197 | CP012599, **CP028197**, CP028537 | 0.842 (CP012599↔CP028537 directly; the group is **chained** through CP028197, whose pairs with each of the other two clear 0.85 at 0.878 and 0.877) |
| 2 | 2 | CP020841 | **CP020841**, CP021960 | 0.995 |
| 3 | 1 | CP006642 | CP006642 | — |
| 4 | 1 | CP010512 | CP010512 | — |
| 5 | 1 | CP017220 | CP017220 | — |
| 6 | 1 | CP018961 | CP018961 | — |
| 7 | 1 | CP023893 | CP023893 | — |
| 8 | 1 | CP026216 | CP026216 | — |
| 9 | 1 | CP026975 | CP026975 | — |
| 10 | 1 | LT991956 | LT991956 | — |
| 11 | 1 | MF344582 | MF344582 | — |

This is not a contradiction of §7 — it is a stricter cross-check. §7 deliberately
restricted the comparison to five conserved genes so locus-private insertions would not
swamp the signal; pangraph's threshold is applied to the whole window, where private
insertion diversity dominates for most loci. The two groups that clear **both** tests —
{CP012599, CP028197, CP028537} and {CP020841, CP021960} — are the strongest candidates
in this pilot for one shared insertion history, having been found independently by two
methods that look at entirely different evidence. The remaining 9 singleton contexts
share no up-/downstream sequence with any other locus at this threshold — including
loci that §7's gene-architecture Jaccard still grouped together (CP023893, CP026975
and LT991956 all landed in §7's "ISEsa1 only" group, but are 3 separate singleton
contexts here): an element catalogue in common without (yet) enough raw sequence in
common to confirm one shared insertion event — exactly the ambiguity §7 warned its own
grouping could not resolve on its own.

An earlier version of this pipeline used a pairwise `blastn -task dc-megablast`
cross-check here instead. It was retired once pangraph's real N-way alignment graph
proved to subsume it — pangraph represents a rearrangement directly, rather than
approximating one by chaining pairwise hits outward from the core.

## 8. The figure

![relBE flanking architecture, aligned on the bait](figures/colour_blocks.png)

Fourteen loci, one row each, aligned on the `relB` start codon (the dashed blue line at
position 0). Reading guide:

- **IS elements**: mustard = forward strand, yellow = reverse strand (which named
  element it is comes from the label above the block, not the colour).
- **Cargo genes**: `relB`/`relE` green/red; `cytb`/`czc` purple/blue (labelled "the
  conserved core"); other cargo in neutral grey.
- **● / —** beside the classification badge marks plasmid vs chromosome.
- **VEHICLE-\* / CORE-\*** two-line badge on the right is the independent two-axis
  classification from §5.
- **Coloured swatch** beside the locus name marks membership in one of the §7
  architecture groups.
- **†** after an IS label marks a target-site duplication shared, letter for letter,
  with another locus (§6).
- The **"idealised architecture" panel** at the bottom is not asserted — every element
  in it is placed at the **median offset** measured across the real loci that carry it,
  and is marked "same in all N" only where that offset is genuinely invariant. `czc` is
  the single fixed point in the whole dataset: **−1184 bp from relB in all nine loci
  that carry it, with zero variation.**

A companion heatmap (`figures/is_matrix.png`) shows IS element presence/absence across
all 14 loci, columns grouped by ISfinder family.

### The same loci, summarised by sequence homology instead

![pangraph: 11 unique contexts summarising 14 loci](figures/pangraph_context_summary.png)

This is §7b's pangraph check, drawn rather than just tabulated: the same 14 loci,
collapsed into their **11 unique up-/downstream contexts**, one row per context, each
drawn from its most-typical ("medoid") member. Reading guide:

- **Orange block at position 0** is the conserved core (`relB`+`relE`), the same anchor
  as the figure above — every context lines up on it.
- **Same colour on either side of the core** (any colour other than pale pink or orange)
  means that exact stretch of DNA is shared with at least one other locus; **different
  colours are different, unrelated combinations of loci** — there is no universal
  meaning to a given hue beyond "look up which loci share it" (see §7b and
  `pangraph_blocks.tsv`'s `group_members` column).
- **Pale pink** is sequence private to that one context — nothing else in the dataset
  shares it.
- **The small swatch beside each context's label** is that context's pangraph group
  (§7b's `pangraph_groups.tsv`) — a different concept from the block colours; it marks
  which contexts cleared the 0.85 shared-sequence threshold together, not which DNA they
  literally share.
- Context 01 (n=3: CP012599, CP028197\*, CP028537) and Context 02 (n=2: CP020841\*,
  CP021960) are exactly the two groups that also independently clustered together in
  §7's gene-architecture Jaccard grouping — the strongest candidates in this pilot for
  one shared insertion history.

Two further figures come from the same pangraph check, one row per
**individual** locus this time rather than per context:

![pangraph: sequence-homology blocks, one row per locus](figures/pangraph_blocks.png)

`pangraph_blocks.png` — same colour scheme as the context figure above (same
colour = mutually alignable, pale pink = private), but every one of the 14
loci gets its own row instead of being collapsed to its context's medoid.

![pangraph background with gene/IS arrows overlaid](figures/pangraph_annotated.png)

`pangraph_annotated.png` — the gene/IS arrows from `colour_blocks.png` (§8,
top) redrawn on top of this sequence-identity background instead of the
flat idealised schematic. The same arrow colour on a *different* background
colour is the same named element sitting on genuinely different DNA — e.g.
`czc` (blue arrow) appears on blue, green, magenta and plain pink
backgrounds across different loci, meaning it is conserved as a *named gene*
far more widely than any single stretch of flanking DNA around it is.

## 9. Next steps: scaling to the whole of NCBI

This pilot deliberately used a small, hand-curated set of 16 records so every step
could be checked by hand against a result verified before the pipeline existed. Scaling
up is not a rewrite — the same Tier A/B/C pipeline and the same annotation layers apply
unchanged — but it does need two things built first, and should proceed in stages
rather than all at once.

1. **Build the extended core anchor (the pivot, not yet done).** The 527 bp `relBE` bait
   alone is too short and too paralogous to anchor a whole-NCBI search: a keyword
   search for `relB`/`relE` returns on the order of 188,000 records, almost all
   unrelated RelE/ParE superfamily background. The fix is to multiple-align the
   `cytb → relE` span across the loci already confirmed in this pilot, find the left
   edge of invariant homology, and write that out as a single extended anchor sequence.
   **Before using it for anything**, that anchor must be BLASTed against generic
   Enterobacteriaceae chromosomal `relBE`/`parDE` homologues to confirm it does **not**
   hit them — if it does, the anchor needs to extend further upstream.
2. **Search by homology, not by keyword**, against a realistic target set — RefSeq
   Enterobacterales complete genomes (~26,500 assemblies) and/or a curated plasmid
   collection (PLSDB) — using the extended anchor from step 1. This step moves to HPC
   (**Gadi or Firefly**, whichever has that database already staged) rather than
   running locally, because a sweep at this scale needs a **pre-built BLAST database
   available on disk**, not fetched per-query over the network — compute-node queues
   on both clusters have no outbound network for anything else, and re-downloading
   ~26,500 assemblies per run is not realistic regardless. Firefly is also where the
   larger annotation tools this pilot could not run locally (ISEScan, IntegronFinder2,
   mob_suite, geNomad) would run at this scale, alongside Bakta.
3. **Stage through the 58-entry TADB3 curated relBE set first**, as a smaller
   intermediate validation step before committing to the full sweep — specifically to
   check what fraction of those 58 carry the full `cytb–czc` core versus a bare TA pair,
   since if most don't, the anchor built in step 1 needs rethinking before it is
   trusted at NCBI scale.
4. **Run the decisive sequence-level test** (all-vs-all `blastn -task dc-megablast`,
   splitting threshold stepped up from the core edge in 10 bp increments, upstream and
   downstream resolved separately, assembled into a `ggtree` splitting tree) once there
   are enough loci for it to carry real statistical weight — the lighter Jaccard
   grouping in §7 is a proxy for this, not a substitute, and explicitly says so.
5. **Cluster by host species vs plasmid backbone** (Jaccard distance + PERMANOVA via R's
   `vegan`) on the resulting IS presence/absence matrix, once the dataset is large
   enough for that test to be meaningful.

Every stage above reuses tooling already built and verified in this pilot — the BLAST
tiering, the fragment-preserving IS layer, the two-axis classification, the gene-family
clustering, and the junction scan all run unchanged on a larger input set. What changes
is only the size of the genome database searched and the anchor used to search it.
