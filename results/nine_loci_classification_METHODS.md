# How the nine Kamaal/Sept-2018 loci were classified

Companion to `results/nine_loci_classification.tsv`.
Run date: 2026-10-26 · BLAST+ 2.17.0 (aarch64) · every command below is reproducible as written.

---

## Inputs

| What | Path | Note |
|---|---|---|
| relBE anchor | `anchors/relBE_seed.fa` | 527 bp; relB 1–240, relE 240–527, 1 bp overlap |
| genome sequences | `data/public_16/fasta/*.fna` | fetched from NCBI, lengths verified vs `LOCUS` |
| genome annotations | `data/public_16/gbk/*.gbk` | all terminate `//` |
| vehicle ends | `/tmp/_vehicle.fa` | `ISEsa1_ISNCY_IS1202` (1609 bp) + `ISEsa2_IS110_unknown` (1493 bp), pulled from `db/isfinder/IS.fna` |

Extracting the two vehicle ends from the ISfinder mirror:

```bash
python3 -c "
s=open('db/isfinder/IS.fna').read().split('>')
keep=[x for x in s if x.split('\n')[0].startswith(('ISEsa1','ISEsa2'))]
open('/tmp/_vehicle.fa','w').write(''.join('>'+k for k in keep))"
```

---

## Step 0 — replicon typing (plasmid vs chromosome)

**Decision rule:** a record is a plasmid **iff** its `source` feature carries a
`/plasmid="..."` qualifier. This is the authoritative GenBank field — not the
DEFINITION line, and not the deck's label.

```bash
grep -m1 -o '/plasmid="[^"]*"' data/public_16/gbk/<ACC>.gbk
```

Present → `plasmid`, and the qualifier value is the plasmid name.
Absent → `chromosome`.

Result: **4 plasmids, 5 chromosomes.** This matches the deck's own split
("RelBE context on plasmids" vs "on chromosomes").

> Naming discrepancy: `MF344582.1` is **pNY2385-HI2** in GenBank, but the deck labels
> it **p525011-HI2**. Same accession, same organism. Worth confirming with Kamaal.
> `CP028197.1` is **pGMI14-002_1**; the deck just says "concord plasmid".

---

## Step 1 — locate relBE

```bash
blastn -task blastn \
       -query anchors/relBE_seed.fa \
       -subject data/public_16/fasta/<ACC>.fna \
       -outfmt "6 pident length sstart send" \
       -evalue 1e-3
```

**Acceptance filter:** `pident >= 90` **AND** `length >= 400`

Rationale: chromosomes carry 5–12 RelE/ParE-family paralogues each (measured). A bare
seed search is not specific. Requiring ≥400 of 527 bp at ≥90% identity demands that
*both* ORFs match, which no superfamily paralogue achieves.

**Strand and anchor point:**

```
  sstart > send   →  strand = minus,  relB_start = max(sstart, send)
  sstart < send   →  strand = plus,   relB_start = min(sstart, send)
```

`relB_start` becomes the **anchor** — offset zero for everything in Step 2.
(Per `CLAUDE.md` this is provisional; final offsets re-zero to the core's left edge
once Stage 2 fixes it.)

All nine passed: eight at **100.0%** identity over the full 527 bp, one
(`CP017220.1`) at **99.4%**.

---

## Step 2 — locate the two vehicle ends

```bash
blastn -task blastn \
       -query /tmp/_vehicle.fa \
       -subject data/public_16/fasta/<ACC>.fna \
       -outfmt "6 qseqid pident length qlen sstart send" \
       -evalue 1e-3
```

`-task blastn` (word size 11), **not** megablast — megablast's word size 28 misses
diverged IS copies, and these sit at 87–92% identity.

Three filters, applied in this order:

| # | Filter | Threshold | Why |
|---|---|---|---|
| 1 | **proximity** | `abs(hit_midpoint - relB_start) <= 40000` | the end must belong to *this* locus. IS elements recur all over a replicon; a copy 2 Mb away is irrelevant. 40 kb is deliberately wider than the ±15 kb analysis window so a distant-but-genuine end is not missed at screening. |
| 2 | **identity** | `pident >= 85` | excludes RelE/IS superfamily background. Genuine ends sit at 87–100%; background at 66–78%. |
| 3 | **length** | `length >= 200` | excludes short spurious alignments. Observed background hits were 64–105 bp at ~70%. |

Where an element survived with several HSPs, the **longest** was kept as representative.

---

## Step 3 — coverage

```
coverage % = 100 * alignment_length / QUERY length
             (qlen = 1609 for ISEsa1, 1493 for ISEsa2)
```

> **Trap, hit for real.** The first run computed `alignment_length / slen`. Because the
> *genome* is the subject, `slen` was the whole replicon (0.3–5.6 Mb), so every element
> reported **0% coverage** and all nine were misclassified. Coverage must be measured
> against the **query** (the IS element), not the subject. If you reuse this code, check
> which way round `-query` and `-subject` are.

---

## Step 4 — the decision rule

```
  ISEsa1 present AND ISEsa2 present
        │
        ├── both coverage >= 90%  ─────────────▶  COMPLETE
        └── either coverage <  90% ────────────▶  DECAYED
  exactly one present ───────────────────────────▶  HALF-CORE
  neither present ───────────────────────────────▶  NO CORE
  (relE-only at <85% id, no relB) ────────────────▶  PARALOGUE  [not triggered here]
```

---

## Results

| # | Accession | Organism | Replicon | relBE | ISEsa1 cov | ISEsa2 cov | Class |
|---|---|---|---|---|---|---|---|
| 1 | CP028197.1 | *S. enterica* Concord | **plasmid** pGMI14-002_1 | 100.0% | 99% | 100% | COMPLETE |
| 2 | CP028537.1 | *E. hormaechei* | **plasmid** pQnrB4_020042 | 100.0% | 99% | 100% | COMPLETE |
| 3 | MF344582.1 | *C. freundii* NY2385 | **plasmid** pNY2385-HI2 | 100.0% | **17%** | 100% | DECAYED |
| 4 | CP018961.1 | *E. coli* Ecol_422 | **plasmid** pEC422_1 | 100.0% | 96% | 100% | COMPLETE\* |
| 5 | CP010512.1 | *E. cloacae* colR/S | chromosome | 100.0% | 99% | **absent** | HALF-CORE |
| 6 | CP026216.1 | *Citrobacter* CFNIH10 | chromosome | 100.0% | 99% | **absent** | HALF-CORE |
| 7 | CP020841.1 | *K. pneumoniae* KPN1482 | chromosome | 100.0% | 99% | 100% | COMPLETE |
| 8 | CP012599.1 | *S.* Newport 0307-213 | chromosome | 100.0% | 99% | 100% | COMPLETE |
| 9 | CP017220.1 | *E. coli* FAM21845 | chromosome | **99.4%** | **absent** | **absent** | NO CORE |

```
  by class      5 COMPLETE · 1 DECAYED · 2 HALF-CORE · 1 NO CORE
  by replicon   4 plasmid  (3 COMPLETE, 1 DECAYED)
                5 chromosome (2 COMPLETE, 2 HALF-CORE, 1 NO CORE)
```

---

## Caveats — read before quoting these classes

1. **"COMPLETE" means all sequence is present, NOT intact or functional.**
   `CP018961.1` is COMPLETE at 96% ISEsa1 coverage, yet that ISEsa1 is **split in two**
   by a 2278 bp ISKpn21 insertion carrying a 5 bp direct repeat. Coverage alone cannot
   see an interruption. A companion `interrupted_by` flag is needed — see PIPELINE.md §4.

2. **"absent" means "not detected under these thresholds", not "deleted".**
   ≥85% identity, ≥200 bp, within 40 kb. A 64 bp ISEsa1 stub was found in pEC422_1,
   so short remnants demonstrably exist and would be invisible here. Before calling
   ISEsa2 genuinely deleted in CP010512/CP026216, run a low-stringency sweep
   (no length cutoff, `-evalue 1e-1`).

3. **The classes are screening calls from two BLAST searches**, not a full structural
   analysis. They do not yet use direct repeats, nesting, or the ±15 kb window.

4. **n = 9.** The plasmid/chromosome pattern below is an observation, not a result.

---

## A pattern worth watching

```
  ISEsa2 present :  4/4 plasmids   vs   2/5 chromosomes
  worst class    :  DECAYED on plasmids   vs   NO CORE on chromosomes
```

The vehicle is better preserved on plasmids than on chromosomes. That is what you would
expect if the plasmid copies are the active, still-mobile form and chromosomal copies
are terminal landings where the transposase ends decay once they are no longer under
selection — consistent with slide 14's proposed direction of travel.

**With n = 9 this is suggestive only** (Fisher's exact on 4/4 vs 2/5 gives p ≈ 0.17).
It becomes testable at the scale of the Enterobacteriaceae expansion.
