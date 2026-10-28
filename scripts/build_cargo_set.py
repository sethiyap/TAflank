#!/usr/bin/env python3
"""
Build a curated cargo-gene nucleotide set from the GenBank records in hand.

  build_cargo_set.py data/public_16/gbk assets/cargo_genes.fa

Why this exists: characterising a locus's context -- which cargo genes sit
beside the bait -- does NOT require a full annotation pass. ISfinder already
covers the IS layer by homology, and a small curated set covers the handful of
cargo genes the biology actually turns on. Both are blastn against something
small, so the whole Tier C classification is cheap.

Bakta is for the figure's complete gene track, and is therefore the LAST and
most expensive step, not a prerequisite for classification.

Sequences are extracted from the project's own records, so every entry is
traceable to an accession and coordinates rather than fetched from elsewhere.
"""

import os, re, sys

# product-name patterns -> the family label used downstream
WANTED = [
    ('czc',       r'cobalt[- ]zinc[- ]cadmium|czc[ABCD]|cation efflux|CzcD|cusA|cusB'),
    ('cytb',      r'cytochrome b/?b6|cytochrome b6|Cytb'),
    ('relB',      r'\brelB\b|antitoxin RelB|RHH.*antitoxin'),
    ('relE',      r'\brelE\b|toxin RelE|mRNA interferase RelE'),
    ('hok_gef',   r'\bhok\b|\bgef\b|Hok/Gef|host killing'),
    ('qnr',       r'\bqnr[A-E]?\d*\b|quinolone resistance'),
    ('tra',       r'\btra[IDJKLMNU]\b|conjugal transfer|relaxase|MobA|TraI'),
    ('intI',      r'integron integrase|\bintI\d?\b'),
    ('mer',       r'\bmer[ARTPC]\b|mercuric'),
    ('ars',       r'\bars[ABCDR]\b|arsenic'),
    ('sil',       r'\bsil[ABCERSP]\b|silver'),
    ('ter',       r'\bter[ABCDEWXYZ]\b|tellurite'),
]
COMP = [(k, re.compile(p, re.I)) for k, p in WANTED]
_T = str.maketrans('ACGTacgtRYKMSWBDHVNrykmswbdhvn',
                   'TGCAtgcaYRMKSWVHDBNyrmkswvhdbn')


def rc(s):
    return s.translate(_T)[::-1]


def parse(path):
    """Yield (product, locus_tag, start, end, is_complement) for each CDS, plus
    the record's sequence and accession."""
    txt = open(path).read()
    acc = re.search(r'^LOCUS\s+(\S+)', txt, re.M).group(1)
    origin = txt.find('\nORIGIN')
    seq = ''.join(re.findall(r'[acgtnACGTN]+', txt[origin:])) if origin > 0 else ''
    feats = txt[:origin] if origin > 0 else txt
    out = []
    for m in re.finditer(r'^     CDS             (\S.*?)(?=^     \S|\Z)',
                         feats, re.M | re.S):
        block = m.group(0)
        loc = m.group(1).split('\n')[0].strip()
        nums = [int(x) for x in re.findall(r'\d+', loc)]
        if not nums:
            continue
        pm = re.search(r'/product="([^"]*)"', block, re.S)
        prod = ' '.join(pm.group(1).split()) if pm else ''
        lt = re.search(r'/locus_tag="([^"]*)"', block)
        out.append((prod, lt.group(1) if lt else '', min(nums), max(nums),
                    'complement' in loc))
    return acc, seq.upper(), out


def main():
    gbkdir, outfa = sys.argv[1], sys.argv[2]
    seen, recs, counts = set(), [], {}
    for f in sorted(os.listdir(gbkdir)):
        if not f.endswith('.gbk'):
            continue
        acc, seq, cds = parse(os.path.join(gbkdir, f))
        if not seq:
            continue
        for prod, lt, lo, hi, comp in cds:
            if not prod:
                continue
            for fam, rx in COMP:
                if not rx.search(prod):
                    continue
                if hi - lo + 1 < 150 or hi > len(seq):
                    continue
                s = seq[lo - 1:hi]
                if comp:
                    s = rc(s)
                key = s[:200]
                if key in seen:          # one representative per distinct allele
                    continue
                seen.add(key)
                counts[fam] = counts.get(fam, 0) + 1
                recs.append((f'{fam}|{acc}|{lt or "na"}|{lo}-{hi}', prod, s))
                break
    os.makedirs(os.path.dirname(outfa) or '.', exist_ok=True)
    with open(outfa, 'w') as fh:
        for name, prod, s in recs:
            fh.write(f'>{name} {prod}\n')
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + '\n')
    print(f'{len(recs)} sequences -> {outfa}')
    for k in sorted(counts):
        print(f'  {k:<10} {counts[k]}')


if __name__ == '__main__':
    main()
