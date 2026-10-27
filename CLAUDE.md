# relBE mobile genetic elements project

Determine structural diversity around relB–relE: use relBE as an anchor to identify
MGEs in the flanking regions, assess conservation of neighbouring genes, and decide
whether relBE arrived in each locus by one ancestral vehicle or by repeated
independent acquisition.

Owner: Pooja Sethiya. Source deck: `relBE_schematic.pptx` (19 slides, dated 19 Aug 2026)
— the scientific rationale lives there; this file is the executable plan and the
record of what has been verified.

Method is modelled on Acman et al. 2022 (blaNDM): anchored, alignment-based
splitting rather than naive flank comparison. Deck slides 5–9 summarise that paper;
slides 11–18 adapt it to relBE.

---

## Project state (verified 2026-10-14)

- This folder contains **only** `relBE_schematic.pptx`. No sequence data, no code yet.
- **Kamaal's nine GenBank sequences are not on this machine** (searched Documents,
  Desktop, Downloads). Slide 11 already shows annotation results for them, so they
  exist elsewhere. They go in `data/kamaal_9/`. Stage 1 onward is blocked without them.
- Tooling on this laptop: only `samtools`, `nextflow`, `R` are installed.
  **Missing:** `blastn`, `makeblastdb`, `seqkit`, `datasets`/EDirect, `bakta`,
  `ISEScan`, `flanker`, `abricate`, `mash`, Biopython.
- `conda` is broken (`__conda_exe: permission denied`). Docker is installed but the
  daemon is not running. Machine is arm64 macOS.
- Open questions for Pooja: is Firefly/HPC access available (needed for the
  big-database tools)? Are raw reads available for any of the nine loci (needed for
  the soft-clip overhang check, Stage 3)?

## The anchor

`anchors/relBE_seed.fa` — the 527 bp sequence Pooja provided. Verified to contain
exactly two ORFs:

| gene | coords (nt) | length | protein |
|---|---|---|---|
| relB (antitoxin, RHH) | 1–240 | 79 aa | `MGSINLRIDDEL…TLDEL` |
| relE (toxin, endoribonuclease) | 240–527 | 95 aa | `MAYFLDFDERAL…VKRIL` |

The ORFs overlap by 1 bp — the classic co-transcribed TA arrangement (slide 10).

**Important:** this 527 bp is the *seed*, not the final anchor. Slide 18 is explicit
that the anchor must be extended through `Cytb/b6 – CO-Zn-Cad (czc) – relB – relE`,
because the TA pair alone is too short and too paralogous to discriminate. The
extended anchor is built in Stage 2 and written to `anchors/relBE_core.fa`.

## What slide 11 already establishes (nine preliminary loci)

- `Cytb/b6` + cobalt–zinc–cadmium resistance genes sit immediately upstream of
  relB–relE in **every** locus. relBE is the tail end of a metal-resistance module.
- **Variation is one-sided**: every insertion is upstream of that core. Only KPN1482
  shows anything downstream at all.
- **ISEsa1-like keeps reappearing**, often degraded — IS903B truncates it in one
  locus, ISKpn38 sits inside it in another. This is the ancestral-vehicle pattern.
- Loci: S. enterica concord plasmid (IncFIB), C. freundii p525011-HI2 (IncFIB),
  E. coli pEC422_1, Citrobacter sp. CFNIH10 (chromosome), E. coli FAM21845
  (chromosome), + others incl. KPN1482.

relBE is **not** itself a mobile element. It is cargo that makes mobile elements
harder to lose (addiction module), which is why it turns up on so many of them.

---

## Repo layout

```
mobile_genetic_elements/
├── CLAUDE.md                 # this file
├── env/relbe.yml             # pinned conda env
├── anchors/
│   ├── relBE_seed.fa         # 527 bp relB+relE
│   └── relBE_core.fa         # extended anchor (Stage 2 output)
├── data/
│   ├── kamaal_9/             # Kamaal's GenBank files, untouched
│   ├── tadb3/                # 58 relBE entries + raw fetch, with fetch date
│   └── genomes/              # downloaded assemblies (gitignore)
├── loci/                     # +/-20 kb extracts, one FASTA per copy
├── results/                  # the tables that are the deliverable
├── scripts/                  # extract_loci.py, relbe_junctions.py, score_mge.py
└── workflow/main.nf          # Nextflow
```

## Conventions

- **Reorient first.** Every window is flipped so relB is on the `+` strand before
  anything else. All offsets depend on this.
- **Coordinate zero is the left edge of the conserved core**, not the relB start codon.
- **One record per copy.** Name loci `{accession}__{contig}__{copy}__{strand}`.
  Paralogous relE-family loci must never collapse into one node (slide 13).
- **+/-20 kb** windows, clamped at contig ends, *recording how much was actually
  obtained on each side*.
- **Never silently drop** contig-truncated loci or fragmentary IS hits. Truncation is
  itself a signal — the blaNDM study lost 2,269 contigs because the insertion hotspot
  breaks assemblies (slide 8). Fragments are the evidence, not noise.
- Follow the `configs/system_paths.yaml` pattern from
  `~/Documents/GitHub/phagpipe/` so paths are not hardcoded.
- Recommended: install **micromamba** rather than repairing anaconda. Run
  `bakta`/`ISEScan`/`geNomad`/`mob_suite` on Firefly (large DBs, x86 wheels); keep
  `blastn`/`seqkit`/R/`gggenes`/`ggtree` local for the iterate-and-look loop.

---

## Stage 0 — environment + scaffold

Create the layout above. Build `env/relbe.yml` with: `blast`, `seqkit`, `samtools`,
`ncbi-datasets-cli`, `entrez-direct`, `biopython`, `pandas`, `bakta`, `isescan`,
`integron_finder`, `mob_suite`, `genomad`, `flanker`, `mash`, `abricate`.

## Stage 1 — one genome, end to end (slide 12)

Use `E. coli pEC422_1` — a closed plasmid, so no contig-edge complication.

1. `makeblastdb`, then `blastn` the 527 bp seed. Confirm one full-length hit; note strand.
2. Reorient so relB is `+`.
3. `samtools faidx` +/-20 kb, clamped, recording actual flank lengths.
4. `bakta` the 40 kb window. Confirm recovery of `Cytb/b6`, `czc`-type CO-Zn-Cad
   genes, then relB/relE.
5. `blastn` the window vs ISfinder. Confirm `ISKpn21` and `ISEcp1-qnrE1` appear, as
   slide 11 claims.

**Gate:** step 5 must reproduce slide 11 for this locus. Fix problems here, not at
500 genomes.

## Stage 2 — build the extended anchor from the nine loci

The pivot of the whole project.

1. Run Stage 1 on all nine. Two are chromosomal (Citrobacter CFNIH10, E. coli
   FAM21845) and will have richer context.
2. Multiple-align the nine `Cytb/b6 → relE` spans. Find the **left edge of invariant
   homology** — that boundary is the top of the conserved core and the origin for
   every offset reported downstream.
3. Write `anchors/relBE_core.fa` (consensus core). Record its length and gene order.
4. **Paralogy check:** `blastn` the core against several RefSeq Enterobacteriaceae
   genomes. Generic chromosomal `relBE`/`parDE` homologues must NOT hit it. If they
   do, extend further upstream.

## Stage 3 — the MGE layer (the actual "find MGEs" step)

Per locus, on the +/-20 kb window, with the core edge as coordinate zero:

| Layer | Tool | Answers |
|---|---|---|
| IS, de novo | `ISEScan` | unbiased IS calls + IR/TSD coordinates |
| IS **fragments** | `blastn` vs ISfinder, **no length cutoff**, report identity *and* query coverage | ISEScan will miss a ~75 bp remnant; the ISEsa1-like stub test needs fragment sensitivity |
| Composite transposons | script: two ISs, same family, same orientation, <30 kb apart, enclosing the core (slide 9 rule) | counts the vehicles; also count how many have *both* flanking ISs intact vs decaying |
| Integrons | `IntegronFinder2` | class 1 integrons are common cargo beside metal-resistance modules — likely relevant to czc |
| Plasmid vs chromosome, plasmid type | `mob_typer` / `mob_recon` | replicon type (IncFIB appears twice on slide 11), relaxase, oriT |
| Prophage | `geNomad` | whether any loci are phage-borne |
| ICE (chromosomal loci) | `ICEfinder` / `oriTfinder` | how relBE entered a chromosome |
| Insertion footprints | `scripts/relbe_junctions.py` — extract +/-25 bp at each MGE boundary, look for direct repeats / TSDs | separates real insertions from assembly artefacts; **distinct DRs across loci = independent events** |
| Insertion footprints, read-level | map raw reads back, screen soft-clipped overhangs >=50 bp against ISfinder (slide 9) | strongest evidence a break is an insertion. **Only possible where raw reads exist** — drop if assemblies only |

**`results/mge_hits.tsv` is the deliverable** (slide 15), one row per MGE call:

```
locus_id, genome, species, replicon, contig, core_edge_pos, element, family,
start, end, strand, offset_from_core, identity, query_cov, is_fragment,
DR_seq, DR_len, truncated_by_contig_edge
```

Companion tables: `loci.tsv` (one row per relBE copy), `core_conservation.tsv`,
`structural_variants.tsv`.

## Stage 4 — TADB3 expansion (58 experimentally validated relBE entries)

Source: https://bioinfo-mml.sjtu.edu.cn/TADB3/browse.php?family=relBE

1. Check for a bulk download link first; scrape `browse.php?family=relBE` only if
   there isn't one. Save raw HTML/TSV to `data/tadb3/` with a fetch date.
2. `efetch`/`datasets` the parent records; extract +/-20 kb per copy using the naming
   convention above.
3. Re-run Stage 3 unchanged. **Record which of the 58 carry the full `Cytb/b6–czc`
   core and which are bare TA pairs** — that split is a result in itself, and shows
   whether the extended anchor generalises or is specific to this metal-resistance
   lineage.

**Gate: review results with Claude/Pooja here before scaling further.** If most of
the 58 lack the core, the anchor strategy needs rethinking before expanding to the
broader GenBank set.

## Stage 5 — structural variants and the decisive test

1. `Flanker` as the quick pass (slide 14) — see where homology drops away.
2. Careful pass: all-vs-all `blastn -task dc-megablast`; keep only hits spanning the
   complete core; cut each hit at the core; **resolve upstream and downstream
   separately**; raise the splitting threshold in 10 bp steps; build the component
   graph; `ggtree` with branch position = distance from core (slide 7).
3. **The decisive test (slide 18):** plot the offset distribution of ISEsa1-like
   fragments relative to the core edge.
   - Fixed offset, present in most loci despite heavy overwriting beyond it →
     **one ancestral vehicle**, everything since is graffiti. (For blaNDM, a ~75 bp
     remnant of ISAba125 in >98% of long contigs identified Tn125 as the ancestral
     vehicle.)
   - No conserved remnant, plus diverse ISs with distinct direct repeats →
     **repeated independent acquisition**.
4. Confirm the one-sided asymmetry from slide 11 holds across all 58, and run the
   downstream analysis anyway so "empty" is demonstrated rather than merely unexamined.

## Stage 6 — grouping by species

Build an IS presence/absence matrix (loci x element). Jaccard distance + PERMANOVA
(`vegan`) with genus/species and replicon type as factors. Answers whether
neighbourhood patterns cluster by host species or by plasmid backbone. Show this
before expanding to the broader GenBank set.

---

## Next action

Pooja is placing Kamaal's nine GenBank files into `data/kamaal_9/`. Once they are
there: scaffold the repo (Stage 0), then run Stage 1 on `E. coli pEC422_1` and check
the gate against slide 11.
