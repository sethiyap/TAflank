# relBE MGE pipeline — filtering and element-calling schematic

Companion to `CLAUDE.md` (which holds the scientific rationale and stage plan).
This file answers two questions only:

1. **How does data get filtered** — what enters, what is kept, what is flagged, what is never dropped.
2. **How are MGEs actually called** — which layer catches which kind of element, and how calls are resolved.

Status key: `[OK]` verified on real data · `[~]` partially verified · `[ ]` planned, not yet run.

---

## Why this project exists

### The observation that started it

From the Sept 2018 lab meeting, slide 1: relBE is a well-studied *E. coli* chromosomal
toxin–antitoxin pair, part of a prophage region with *hokD* (formerly *relF*) downstream.
It also turns up on plasmids, nearly identical in sequence. Then:

> *"Kamal noticed apparently related transposase genes annotated around relBE in
> different locations. **Mobilisation of relBE?**"*

That is the whole project. Transposases keep appearing beside relBE in unrelated
organisms and on unrelated replicons. Either something is carrying relBE around, or
relBE keeps being picked up independently. We are here to tell which.

### Why relBE is worth chasing

relBE is a **type II toxin–antitoxin addiction module**, and that makes it a very
particular kind of passenger:

```
   RelE  =  toxin,     stable,   mRNA interferase (cleaves mRNA in the ribosomal A site)
   RelB  =  antitoxin, LABILE,   neutralises RelE — but is constantly degraded

   cell keeps the replicon  →  RelB replenished  →  RelE neutralised  →  cell lives
   cell LOSES the replicon  →  RelB decays       →  RelE free        →  cell dies/arrests
```

This is **post-segregational killing**. A replicon carrying relBE is punishingly hard to
lose: daughter cells that drop it are removed from the population.

So the key framing — and the reason the analysis is shaped the way it is:

> **relBE is not itself a mobile element. It is cargo that makes mobile elements harder
> to lose.** That is precisely why it turns up on so many of them.

We are not tracking a selfish element. We are tracking a **stabilisation function**
being delivered onto things that benefit from persisting.

### Why it matters beyond this gene

Look at what relBE is attached to in every one of our sixteen records:

| Neighbour | What it is | Why it matters |
|---|---|---|
| `czc` / cation-efflux | cobalt–zinc–cadmium resistance | **metal** resistance |
| `qnrE1` | fluoroquinolone resistance | antibiotic resistance |
| IncHI2 backbone | classic MDR plasmid lineage | carries many AMR genes |
| pNDM-CIT-type | carbapenemase-associated platform | last-line resistance |

relBE sits at the tail end of a **metal-resistance module** on **multidrug-resistance
plasmids**. If relBE stabilises these replicons, it contributes to a well-known and
troubling pattern:

```
   metals in the environment          (agriculture, aquaculture, sewage, industry)
            │   no antibiotics needed for this step
            ▼
   selection for czc  ─────────────▶  module retained
            │
            ├── relBE rides along on the same module
            │        └──▶ post-segregational killing → replicon cannot be lost
            │
            └── AMR genes on the same plasmid persist as a by-product
                         ▲
                 CO-SELECTION: resistance maintained with no antibiotic pressure at all
```

That is the public-health relevance: it helps explain why MDR plasmids stubbornly
persist in populations even when antibiotic use is reduced.

### What we are actually looking at — three nested levels

Keeping these separate is what stops the analysis becoming muddled, because the answer
is different at each level:

```
 ┌─ LEVEL 3 ── HOST REPLICON ──────────────────────────────────────────────┐
 │   IncHI2 / IncFIB plasmid, or a chromosome                              │
 │   Q: did the whole replicon move, or just a piece of it?                │
 │                                                                          │
 │   ┌─ LEVEL 2 ── THE VEHICLE ───────────────────────────────────────┐    │
 │   │   ISEsa1b ── ISEsa2 ── [czc cargo] ── relBE,  in a ~16 bp DR    │    │
 │   │   Q: one ancestral transposition unit, or many?                 │    │
 │   │                                                                  │    │
 │   │   ┌─ LEVEL 1 ── THE CARGO ──────────────────────────────┐       │    │
 │   │   │   relB + relE, 527 bp, 1 bp ORF overlap             │       │    │
 │   │   │   Q: one lineage, or convergent recruitment of      │       │    │
 │   │   │      different relBE paralogues?                    │       │    │
 │   │   └──────────────────────────────────────────────────────┘       │    │
 │   └──────────────────────────────────────────────────────────────────┘    │
 └──────────────────────────────────────────────────────────────────────────┘
```

Level 1 is the anchor. Level 2 is where the DR evidence decides the question.
Level 3 is where host range and the species/backbone grouping live.

### A directional clue we have already seen

Slide 1 records that in the *E. coli* chromosome, relBE sits in a prophage region with
*hokD*/*relF* immediately downstream. In pEC422_1 — a **plasmid** — we found a
**Hok/Gef family protein at +72 bp downstream of relE**, exactly that arrangement
preserved.

If the plasmid copies retain the ancestral chromosomal neighbour, that is evidence for
slide 14's proposed direction of travel:

```
   E. coli chromosome  ──▶  plasmids  ──▶  other species' chromosomes
```

rather than the reverse. Testing whether *hok/gif* is retained across all loci is a
cheap, concrete check on that hypothesis.

---

## 0. What we are looking for

The target is not relBE. The target is the **vehicle** that carries it.

From the Sept 2018 deck (spatial layout of slides 7/9/11), the unit is two-ended and
bracketed by a long direct repeat:

```
 [DR ~16bp]                                                                 [DR ~16bp]
     │                                                                           │
     ├── ISEsa1b ──── ISEsa2 ──┬─ Cytb/b6? ─┬─ czc ─┬─ relB─relE ─┬─ hok/gif ─────┤
     │   ISNCY        IS110    │ (not in    │cation-│  527 bp seed│               │
     │   IS1202       no IRs   │  all loci) │efflux │             │               │
     │                         │            │       │             │               │
     └─── VEHICLE ─────────────┴──────────────── CARGO ───────────────────────────┘
          encodes mobility      └─ PAYLOAD ─────┘   └─── ADDICTION ──────┘
          (transposases)          acts on the HOST     acts on the REPLICON
```

**Everything that is not a transposase is cargo** — `Cytb/b6`, `czc`, relBE and
`hok/gif` alike. Only ISEsa1b and ISEsa2 encode mobility; they are the vehicle. But
cargo is not one thing, and the split matters:

| class | examples | acts on | why it persists |
|---|---|---|---|
| **payload** | `czc`, `Cytb/b6`, `qnrE1` | the **host cell** — confers a phenotype | selection favours it directly: metal or antibiotic survival |
| **addiction** | **relBE**, `hok/gif` | the **replicon** — confers nothing on the host | post-segregational killing: cells that lose it die |

So `czc` explains why the module is selected **for**; relBE explains why it cannot be
selected **against**. Two different mechanisms, both riding the same vehicle.

**They are separable — this is measured, not assumed.** Across `data/public_16/`:

```
  czc present, relBE absent  :  5 records  (CP006636, CP018962, CP023893,
                                            CP026975, MF344582)
  relBE present, czc absent  :  1 record   (CP028197)
```

So this is not one indivisible module. It is a vehicle that has accumulated both kinds
of cargo, and the two must be tracked **independently**.

> **Pipeline consequence.** The extended anchor must **not** require `czc`. Defining the
> core as "czc + relBE" would have discarded CP028197 outright. `czc` is context to be
> recorded per locus, not a membership criterion. Same for `Cytb/b6`, which is absent
> from pEC422_1 altogether (see F3).

Two consequences that drive every filter below:

- **The IS ends are annotated under generic family names.** GenBank calls ISEsa1b an
  "ISNCY family transposase" and ISEsa2 an "IS110 family transposase". ISEcp1 appears
  only as "IS1380 family transposase". Annotation keywords cannot be trusted to find
  the vehicle — only sequence homology can.
- **relBE sits at the downstream edge of the vehicle.** Most structural variation is
  therefore upstream. "Most", not "all" — see F6.

### So what is the MGE here?

The relBE **operon itself is not an MGE**. It has no transposase, no inverted repeats,
no recombinase, no origin of transfer. It cannot move itself. Neither can `czc`.

The MGEs are the things *around* it, and they sit at four nested levels:

```
  LEVEL D   chromosomal integration
            whole/partial plasmid lands in a chromosome
            evidence: slide 13 — CP012599, CP010512, CP026975
    │
    └─ LEVEL C   the plasmid              ← an MGE: moves between CELLS
         IncHI2 / IncFIB, conjugative
         evidence: traI, traD, traJ found in the pEC422_1 window
         │
         └─ LEVEL B   the transposition unit (TPU)   ← THE TARGET
              ISEsa1b ─ ISEsa2 ─ [cargo incl. relBE], in a ~16 bp DR
              moves between REPLICONS as one piece
              │
              └─ LEVEL A   individual IS elements   ← MGEs: self-mobile
                   ISEsa1b, ISEsa2, ISKpn21, ISKpn28, ISEcp1,
                   ISKpn38, IS903, ISCfr4, ISSen4, ISKpn8
                   each encodes its own transposase
```

**We are looking for Level B.** Levels A, C and D are observable with standard tools.
Level B is the one that has to be *reconstructed*, because nothing annotates it — and it
is the level that answers the project's question.

> **A precision worth keeping.** Level B is **not a classic composite transposon**.
> In Tn5 or Tn10 the cargo is flanked by two copies of the *same* IS. Here the two ends
> are *different* families — ISEsa1b (ISNCY/IS1202) and ISEsa2 (IS110). The deck's own
> term is the right one: a **transposition unit (TPU)**, as in slide 13's
> "ISEcp1-qnrE1 TPU" and slide 14's "relBE mobilised by ISEsa1-like element in TPU?".
> ISEcp1-style one-ended transposition can mobilise adjacent DNA using one IS end plus
> a fortuitous end in the neighbour, which produces exactly this architecture and
> explains slide 14's "insertion made up of 2 components".
>
> This matters for Layer 3: a composite-transposon rule requiring *same family, same
> orientation* will **not** fire on this structure. It needs a two-different-ends rule.

### What we are looking for, in one line

> Is there **one** TPU — with consistent ends and a shared ~16 bp DR — that delivered
> relBE into this lineage, or did relBE get picked up **many times** by different
> elements leaving different DRs?

---

## 1. Pipeline overview

### 1a. The short version — which tool runs where

| # | Question asked | Tool | In → Out | Status |
|---|---|---|---|---|
| 1 | Where is relBE in this genome? | **`blastn`** (seed vs genome) | genome → coords + strand | `[OK]` |
| 2 | Cut out the neighbourhood | **`samtools faidx`** + revcomp | genome + coords → 30 kb window | `[OK]` |
| 3 | What genes are in the window? | **GenBank annotation**; `bakta` where absent | window → gene list | `[~]` |
| 4 | Which **named IS** elements are here? | **`blastn` vs ISfinder**, no length cutoff | window → IS hits + coverage | `[OK]` |
| 5 | Any IS the database missed? | **`ISEScan`** (de novo) | window → IS + IR/TSD | `[ ]` |
| 6 | Do two IS bracket the cargo? | **custom script** (two-different-ends rule) | IS hits → TPU call | `[ ]` |
| 7 | Is there an integron? | **`IntegronFinder2`** | window → integron | `[ ]` |
| 8 | Plasmid or chromosome, which type? | **`mob_typer` / `mob_recon`** | replicon → Inc type, oriT | `[ ]` |
| 9 | Is it phage-borne? | **`geNomad`** | window → prophage | `[ ]` |
| 10 | How did it enter a chromosome? | **`ICEfinder` / `oriTfinder`** | chr window → ICE | `[ ]` |
| 11 | **Is there a direct repeat at the junction?** | **custom script** (mismatch-tolerant) | window + IS coords → DR | `[OK]` |
| 12 | How do the loci relate to each other? | **`blastn` all-vs-all** `dc-megablast` + R | all loci → graph, `ggtree` | `[ ]` |
| 13 | Do they cluster by host or by backbone? | **R `vegan`** (Jaccard + PERMANOVA) | IS matrix → stats | `[ ]` |

**BLAST does four distinct jobs** — worth keeping separate in your head:

```
  blastn #1   seed vs genome          "where is relBE?"            → step 1
  blastn #2   window vs ISfinder      "which IS are here?"         → step 4
  blastn #3   core vs RefSeq          "is my anchor specific?"     → Stage 2 paralogy check
  blastn #4   loci vs loci (all-vs-all) "how do loci relate?"      → step 12
```

Steps 1, 2, 4 and 11 are the **core loop** and all work today. Everything else is either
an extra evidence layer (5–10) or downstream synthesis (12–13).

### 1b. The filter flow in detail

```
  ┌────────────────────────────────────────────────────────────────┐
  │ INPUT                                                          │
  │   data/public_16/   16 closed GenBank records (Sept 2018 deck) │
  │   data/tadb3/       58 TADB3 relBE entries        [ ]          │
  │   data/genomes/     broader GenBank sweep         [ ]          │
  └───────────────────────────┬────────────────────────────────────┘
                              ▼
       ╔══════════════════════════════════════════════╗
       ║  F1  ANCHOR SEARCH    seed → candidate copies ║   filters out paralogues
       ╚══════════════════════╤═══════════════════════╝
                              ▼
       ╔══════════════════════════════════════════════╗
       ║  F2  REORIENT + WINDOW   ±15 kb, relB on (+)  ║   no filtering; records truncation
       ╚══════════════════════╤═══════════════════════╝
                              ▼
       ╔══════════════════════════════════════════════╗
       ║  F3  CORE TEST        vehicle present or not? ║   splits core+ / bare-TA
       ╚══════════════════════╤═══════════════════════╝
                              ▼
       ╔══════════════════════════════════════════════╗
       ║  F4  MGE CALLING      layered, see §3         ║   never drops fragments
       ╚══════════════════════╤═══════════════════════╝
                              ▼
       ╔══════════════════════════════════════════════╗
       ║  F5  TILING           resolve overlapping hits║   detects nesting
       ╚══════════════════════╤═══════════════════════╝
                              ▼
       ╔══════════════════════════════════════════════╗
       ║  F6  JUNCTIONS        direct repeats / TSDs   ║   the decisive evidence
       ╚══════════════════════╤═══════════════════════╝
                              ▼
                    results/*.tsv  →  §5 decision
```

---

## 2. The filter cascade

### F1 — Anchor search `[OK]`

```
query: anchors/relBE_seed.fa (527 bp)      blastn -evalue 1e-5
  │
  ├── pident ≥ 95, qcov ≥ 95 ────────────────▶ FULL relBE       → carry forward
  │
  ├── relB half only, pident ≥ 95 ───────────▶ ORPHAN ANTITOXIN → carry + flag
  │      harmless, possibly adaptive (immunity to incoming TA systems)
  │
  ├── relE half only ─── pident ≥ 95 ────────▶ ORPHAN TOXIN     → carry + FLAG LOUDLY
  │      │                  our lineage, lost relB — see ORF integrity check
  │      └─── pident < 85 ───────────────────▶ PARALOGUE        → log, do not carry
  │                           RelE/ParE superfamily background
  │
  └── below thresholds ──────────────────────▶ reject (logged, never silently)
```

**Why the relE-only branch splits by identity.** A relE-only hit is one of two completely
different things, and identity separates them cleanly: ~99% means a decayed member of
*our* lineage that lost its antitoxin; 70–85% is RelE/ParE superfamily background, of
which chromosomes carry 5–12 copies each.

### ORF integrity check — mandatory for any orphan toxin

An expressed, functional relE with no antitoxin arrests translation: the cell stops
growing, forms no colony and is never sequenced. **So an orphan toxin in a closed genome
must be explained**, and the check is cheap:

```
for an orphan relE, translate and report:
    internal stop codon?      → pseudogene
    frameshift vs reference?  → pseudogene
    5' or 3' truncation?      → fragment
    intact ORF, full length?  → *** investigate individually ***
                                 promoter lost with relB (silent)?
                                 cross-neutralised in trans by a paralogue?
                                 or mis-annotation?
columns: orf_intact, internal_stop, frameshift, truncated_end, pct_len_vs_ref
```

> **Asymmetry worth remembering.** Orphan **relB** is biologically unremarkable — an
> antitoxin with no toxin is inert, and can even protect against incoming TA systems.
> Orphan **relE** should be close to unobservable. The 1 bp `ATGA` overlap means relE
> translation is coupled to relB's, so simple mutation rarely produces one.
>
> Checked in `data/public_16/`: every genuine RelE has an antitoxin immediately
> adjacent — annotated as "bifunctional antitoxin/transcriptional repressor" or
> "DNA-binding protein", **never as `relB`** — usually with the same 1 bp overlap.
> **No true orphan toxin in our set.** Finding one would be the interesting case.

**Why the paralogue branch matters.** Grepping the 16 GenBank records showed chromosomes
carry **5–12** annotated relE-family genes each (CP026216 has 12, CP020841 has 11). The
bare seed is *not* specific. Full-length 95/95 hits are the discriminator; `relE`-only
hits are RelE/ParE superfamily noise.

> Per `CLAUDE.md`: paralogous relE-family loci must never collapse into one node.
> Each surviving copy gets its own record: `{accession}__{contig}__{copy}__{strand}`.

Verified on CP018961: exactly **1** hit, 100% id, 527/527, minus strand. No ambiguity.

### F2 — Reorient and window `[OK]`

```
if sstart > send:  strand = minus  →  reverse-complement the window
                                      (all downstream offsets assume relB on +)

window = [relB_start − 15000, relE_end + 15000]   clamped to contig ends
record: flank_up_actual, flank_dn_actual, truncated_by_contig_edge (yes/no/left/right)
```

**Why ±15 kb and not ±10 kb.** Measured, not guessed. Offsets of every genuine
ISEsa1/ISEsa2 hit (≥85% id, ≥200 bp) from relB start across the 16 records:

```
  furthest genuine vehicle edge observed:  -12,244 bp   (CP023893, CP026975, LT991956)

   +/- 5 kb : 11/19 vehicle ends CLIPPED
   +/-10 kb :  4/19 vehicle ends CLIPPED   <-- loses the left end in 4 loci
   +/-15 kb :  0/19 clipped                <-- chosen: safe, +2.8 kb margin
   +/-20 kb :  0/19 clipped                    (old default, 25% wasted)
```

±10 kb would truncate the vehicle's ISEsa1 end in CP018961 (-10554) and in the whole
pECL_A group (-12244). ±15 kb keeps every observed end with headroom, at 25% less
sequence than the old ±20 kb.

This step **filters nothing**. Truncation is recorded, never used to exclude — the
blaNDM study lost 2,269 contigs exactly this way. A locus with 3 kb of left flank is
still a locus; it just answers fewer questions.

Verified on CP018961: minus strand → revcomp → relBE lands at exactly 20001–20527,
15000 bp both sides, no clamping (circular replicon).

### F3 — Core test `[~]`

### What "the core" is

```
        <────────────────────── THE CORE ──────────────────────>

   [DR]   ISEsa1b   ISEsa2   [ Cytb/b6? | czc ]   relB─relE   hok/gif   [DR]
          ~1.6 kb   ~1.5 kb      optional cargo     527 bp     optional
          │         │                               │
          └ left    └ right end                     └ ANCHOR (offset zero
            end       of vehicle                      is measured from here)

   canonical offsets from relB start, measured across the 16 records:
          ISEsa1  -5494 .. -3899        ISEsa2  -3898 .. -2406
                                                          ▲
                             ISEsa2 ends at -2406 in SIX independent loci
                             (CP012599, CP020841, CP021960, CP028197,
                              CP028537, MF344582) — a FIXED OFFSET
```

Minimum definition: **ISEsa1b + ISEsa2 upstream of an intact relBE.** `czc` and
`Cytb/b6` are recorded as context, **not** required for membership (see §0).

### The classification — every relBE copy lands in exactly one class

| class | definition | what it means | observed in `data/public_16/` |
|---|---|---|---|
| **COMPLETE** | ISEsa1 **+** ISEsa2, both ≥90% coverage, intact relBE | vehicle intact; the reference architecture | CP012599, CP020841, CP021960, CP028197, CP028537 |
| **DECAYED** | both ends present but one truncated or interrupted | vehicle can no longer move; cargo stuck in place — *this is the informative class* | CP018961 (ISKpn21 inserted into ISEsa1), MF344582 (ISEsa1 only 275 bp left) |
| **HALF-CORE** | exactly one vehicle end present | either a different architecture, or the other end fully lost | CP010512, CP026216, CP023893, CP026975, LT991956 |
| **NO CORE** | relBE present, **no** vehicle end at all | bare TA pair — relBE arrived some other way, or the vehicle is gone | CP006642, CP017220 |
| **PARALOGUE** | relE-family homologue, <85% id, not our lineage | superfamily background — logged, not carried | 56/57 TADB3 loci |

```
   COMPLETE    [DR] ISEsa1 ─ ISEsa2 ─── czc ─── relBE        vehicle intact
   DECAYED     [DR] ISEsa1╪IS╪ ─ ISEsa2 ─ czc ─ relBE        ends broken, cargo stays
   HALF-CORE   [DR] ISEsa1 ──────────────────── relBE        one end only
   NO CORE                                      relBE        bare TA pair
   PARALOGUE                                    relE'        different lineage entirely
```

**Architecture sub-classes** (by ISEsa1 offset — these look like genuinely different
insertion events, not noise):

```
  Class I   canonical HI2   ISEsa1 -5494..-3899 + ISEsa2 -3898..-2406   6 loci
  Class II  pECL_A group    ISEsa1 -12244..-10639, no ISEsa2            3 loci
            (CP023893, CP026975, LT991956 — identical offsets)
  Class III other           CP010512 -5224..-3633 ; CP026216 -3991..-2393
```

The COMPLETE/DECAYED/HALF/NO-CORE split **is itself a result** — it says whether the
vehicle generalises or is specific to this metal-resistance lineage. Across TADB3 the
answer was emphatic: 56/57 are PARALOGUE, 1 is our lineage.

> Open issue: `CLAUDE.md` asserts `Cytb/b6` is present upstream in *every* locus.
> It is **absent from pEC422_1** — no cytochrome b/b6 anywhere in the 30 kb window,
> and slide 5 only draws it on three of the nine. The core definition must be settled
> by homology across all nine before this filter can be finalised.

---

## 3. How MGEs are found (F4)

Layered, because no single tool sees everything. Each layer has a known blind spot.

```
 LAYER 1  de novo IS           ISEScan              IR/TSD coords, unbiased
          blind to: elements with no inverted repeats, and short remnants
          ⚠ ISEsa2 has NO IRs → ISEScan will miss the vehicle's right end

 LAYER 2  IS by homology       blastn vs ISfinder   NO LENGTH CUTOFF        [OK]
          -task blastn (word size 11) for diverged copies
          -evalue 1e-3, -max_target_seqs 5000
          catches: the 64 bp ISEsa1 remnant that Layer 1 cannot see

 LAYER 3  composite Tn         script               two IS, same family,
          same orientation, <30 kb apart, enclosing the core
          also: how many retain BOTH flanking IS intact vs decaying

 LAYER 4  integrons            IntegronFinder2      class 1 beside czc?
 LAYER 5  plasmid context      mob_typer/mob_recon  replicon, relaxase, oriT
 LAYER 6  prophage             geNomad              any phage-borne loci?
 LAYER 7  ICE (chromosomal)    ICEfinder/oriTfinder how relBE entered a chromosome
 LAYER 8  junctions            scripts/relbe_junctions.py   §4
 LAYER 9  read-level           soft-clip ≥50 bp vs ISfinder — only where reads exist
```

**Classification of every Layer-2 hit by subject coverage:**

| coverage | class | treatment |
|---|---|---|
| ≥ 90% | `complete` | named element, full confidence |
| 50–90% | `partial` | named, flagged — likely truncated by a later event |
| < 50% | `fragment` | **kept** — these are the ancestral-vehicle evidence |

Fragments are never filtered out. The whole Stage 5 test depends on a ~64 bp remnant.

### Database note

ISfinder forbids bulk download. `db/isfinder/IS.fna` (5,970 sequences) is the community
mirror — see `db/isfinder/PROVENANCE.txt`. All ten elements named in the 2018 deck are
present. **Re-verify name-level calls on the ISfinder web BLAST before publishing.**

---

## 4. Tiling and nesting (F5) — and one trap

Raw blastn returns overlapping HSPs. Collapsing them needs care.

```
sort HSPs by bitscore descending
for each HSP:
    overlap with any kept interval > TOLERANCE ?  → discard
    else                                         → keep
TOLERANCE = 30 bp
```

> **Trap, hit for real.** With `TOLERANCE = 0`, a **3 bp** overlap discarded a **1540 bp**
> ISEsa1 alignment, making a near-complete element look like a bare 64 bp stub and
> inverting the biological reading. Adjacent elements legitimately share a few bases at
> TSDs and junctions. Never use zero tolerance.

### Nesting detection

```
IF element A yields two HSPs split in the window
   AND their SUBJECT ranges are near-contiguous (A bases 1-1540 and 1533-1596)
   AND element B lies between them in WINDOW coordinates
   AND the shared subject range appears twice in the genome
THEN  B is nested in A; the shared range is the target site duplication
      → mge_hits.tsv: nested_in = A, interrupts = A, DR_seq, DR_len
```

Verified on CP018961: ISEsa1 split into bases 1–1540 and 1533–1596, with ISKpn21
(2278 bp, 100% id) between them, sharing 8 bp — exactly slide 13's
*"ISKpn21 inserted in ISEsa1-like with DR"*.

---

## 5. Direct repeats (F6) — the decisive measurement

Exact-match scanning is **not** sufficient: a TSD that acquired one point mutation since
insertion becomes invisible. Scan with mismatches and read the identity profile.

```
for k in 4..20:
    for boundary shifts −10..+10 on each side of the element:
        compare left-flank k-mer vs right-flank k-mer
        keep best identity
report the full k vs %identity profile; the DR length is where identity COLLAPSES
```

Observed on CP018961 at the ISKpn21 junction:

```
  k :  4     5     6     7     8  │  9    10    11    12
 %id: 100   100   83    86   87.5 │ 78    70    64    58
                            ▲     │  ▲
                   real DR ─┘     └── collapse → boundary is here
```

5 bp exact (`TGAGC`), 8 bp with one mismatch (`TGAGCNGT`). The 8 bp reading agrees
independently with the blast-derived shared subject range — two methods, same answer.

### What DR length means

| DR | interpretation |
|---|---|
| **~16 bp, `AGGA…` family, shared across loci** | the ancestral vehicle. 16 bp cannot recur by chance. |
| 4–9 bp, element-specific, differing between loci | later independent IS insertions — graffiti |
| none | assembly artefact, or a non-TSD-generating mechanism — flag, do not assume insertion |

---

## 6. The decision (Stage 5)

```
                 Is the ~16 bp AGGA-family DR shared across loci?
                        │                           │
                       YES                          NO
                        │                           │
          ONE ANCESTRAL VEHICLE          Are DRs diverse + element-specific?
          everything later is graffiti            │            │
                        │                        YES           NO
                        ▼                         │            │
          Do chromosomal loci differ      REPEATED INDEPENDENT  inconclusive →
          in HOW they arrived?            ACQUISITION           extend anchor,
                        │                                       re-run
                       YES
                        ▼
        ╔══════════════════════════════════════════════════════╗
        ║  BOTH, at two timescales:                            ║
        ║   · one vehicle assembled relBE into the module      ║
        ║   · that module then moved into chromosomes          ║
        ║     repeatedly, by ≥3 different routes (slide 13)    ║
        ╚══════════════════════════════════════════════════════╝
```

Slide 13 already assigns three distinct chromosomal entry routes — direct ISEsa1b
(CP026216), ISEcp1-mediated (CP020841/CP021960), and whole-plasmid integration
(CP012599/CP010512/CP026975). The question is not either/or.

**Also run the downstream side explicitly.** `CLAUDE.md` says only KPN1482 shows anything
downstream, but pEC422_1 carries a second complete ISEcp1 (100% id, 100% cov) at
**+6028, downstream of relBE**. "Empty" must be demonstrated, not assumed.

---

## 7. Output schemas

**`results/loci.tsv`** — one row per relBE copy
```
locus_id, genome, species, replicon, replicon_type, contig, contig_len, topology,
relB_start_orig, relE_end_orig, strand_orig, reoriented, window_start, window_end,
window_len, flank_up, flank_dn, truncated, relBE_in_window,
seed_pident, seed_qcov, n_seed_hits
```

**`results/mge_hits.tsv`** — one row per MGE call *(extended beyond CLAUDE.md with the
last five columns, because the evidence is nested, not flat)*
```
locus_id, element, family, win_start, win_end, len, identity, query_cov, is_fragment,
offset_from_relB, nested_in, interrupts, DR_seq, DR_len, note
```

Companion tables: `core_conservation.tsv`, `structural_variants.tsv`.

> Offsets are currently relative to **relB start**. After Stage 2 fixes the core's left
> edge, every offset is re-zeroed to that edge per the `CLAUDE.md` convention.

---

## 8. Current state

| Stage | Status |
|---|---|
| 0 — environment, scaffold | `[~]` BLAST+ 2.17.0 arm64 at `~/opt/`; ISfinder mirror local. bakta/ISEScan/geNomad/mob_suite still absent |
| 1 — one genome end to end | `[OK]` **gate passed** on pEC422_1 (CP018961) |
| 2 — extended anchor from nine loci | `[ ]` needs Stage 1 on remaining eight |
| 3 — MGE layer | `[~]` Layer 2 working; Layers 1, 3–9 pending |
| 4 — TADB3 expansion | `[ ]` |
| 5 — structural variants, decisive test | `[ ]` |
| 6 — grouping by species | `[ ]` two tests ready: 0/37 vs 51/108 E. coli HI2 split (Fisher); chromosome co-carriage |

**Data in hand:** 16 closed GenBank records, all lengths verified against their `LOCUS`
lines, all terminating `//`. Two negative controls confirmed carrying zero `relB`
(CP018962, CP006636) while their own plasmids carry one.

**Known gap:** `LT991956.1` arrived from ENA with **0 CDS features** — needs annotation
before Layer 3+ can run on it.

**Blocked on tooling:** Layers 1 and 3–7 need bakta / ISEScan / IntegronFinder2 /
mob_suite / geNomad. conda is broken; micromamba or Firefly is the route (see `CLAUDE.md`).


---

## 9. What success looks like

### The four concrete outputs

**1. A delimited, named element — or a demonstration that there isn't one.**
If the ~16 bp `AGGA…` DR brackets the same ISEsa1b–ISEsa2–czc–relBE span across loci,
we can define that element: its exact ends, its DR, its cargo, its length. That is a
nameable unit, the way **Tn125** became the nameable vehicle for *bla*NDM. If it does
not hold up, the negative result is equally publishable and redirects the question.

**2. The tables** (`loci.tsv`, `mge_hits.tsv`, `core_conservation.tsv`,
`structural_variants.tsv`) — a reusable census of what sits around relBE, with every
fragment and truncation retained rather than filtered away.

**3. The figure set:**

| Figure | Shows | Status |
|---|---|---|
| Gene-arrow maps, all loci, aligned on the core | the one-sided/two-sided structure | `gggenes` |
| ISEsa1 fragment offset distribution vs core edge | fixed offset ⇒ ancestral vehicle | the slide-18 test |
| DR alignment across loci | 16 bp conserved vs element-specific short DRs | the decisive panel |
| IS presence/absence matrix + ordination | clustering by species vs by backbone | Jaccard + PERMANOVA |

**4. Testable predictions for the bench.** Slide 14 ends on *"Have at least 3 local
plasmids with relBE. Experiments?"* The analysis should hand the wet lab specific,
falsifiable claims — see below.

### What each possible answer would mean

| Finding | Interpretation | So what |
|---|---|---|
| 16 bp DR shared across loci, ISEsa1 remnant at a fixed offset | **one ancestral vehicle**; all later IS are graffiti | the module is a single traceable unit; it can be named, dated and used as a surveillance marker |
| DRs diverse and element-specific; no conserved remnant | **repeated independent acquisition** | relBE is being recruited again and again ⇒ strong ongoing selection for plasmid stabilisation in this lineage; watch the *function*, not one element |
| Both, at different levels *(currently expected)* | one vehicle assembled czc+relBE; that module then entered chromosomes **repeatedly**, by ≥3 routes (slide 13) | explains why relBE looks both highly conserved and widely scattered |
| Most TADB3 entries are bare TA pairs without the core | the vehicle is specific to this metal-resistance lineage | narrows the claim, and is a result in itself (Stage 4 gate) |

### Predictions the bench could test

If the one-vehicle model holds, these follow and are checkable:

1. **Deleting relBE raises plasmid loss rate.** Stability assay on the three local
   relBE plasmids, ±relBE, passaged without selection. The model predicts measurable
   segregational loss only in the deletion.
2. **Cadmium/zinc exposure maintains the AMR plasmid without antibiotics.** Direct test
   of the co-selection chain above.
3. **The 16 bp DR marks a target site.** If insertion is targeted rather than random,
   the same motif should be findable at unoccupied sites in relBE-free relatives.
4. **Host-range prediction.** `0/37` IncHI2 plasmids with relBE are from *E. coli*,
   versus `51/108` without (slide 7). If real, *E. coli* should resist acquiring or
   retaining relBE-bearing HI2 plasmids — testable by conjugation into *E. coli*.
   This is a Fisher's exact test we can run today on existing data.

### Success criteria for the pipeline itself

The pipeline has done its job when all of these are true:

- [ ] Every relBE copy has one record, with paralogues kept separate, never collapsed
- [ ] Every IS fragment is retained with identity **and** coverage — no length cutoff anywhere
- [ ] Every truncated locus is flagged, not dropped; "empty" regions are demonstrated, not assumed
- [ ] Core edge fixed, and every offset re-zeroed to it
- [ ] DR called for each junction with its identity profile, not a bare sequence
- [ ] Nesting recorded (`nested_in` / `interrupts`), so element-in-element history is readable
- [ ] The §6 decision tree resolves to one branch, with the evidence that put it there
- [ ] Name-level IS calls re-verified on ISfinder web BLAST before publication

### The one-sentence version

> **relBE is a molecular seatbelt.** We are working out whether it was fitted once, by a
> single transposition unit that we can name and delimit, and then inherited — or whether
> it keeps being bolted on independently because anything carrying metal- and
> antibiotic-resistance genes is under strong selection not to be lost.
