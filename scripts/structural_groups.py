#!/usr/bin/env python3
"""
Stage 6 (first half) -- group loci by how similar their core architecture is,
as a direct answer to "how many of these plausibly share one history".

  structural_groups.py --blocks run_all/blocks --out results/structural_groups.tsv

WHY THE VOCABULARY IS RESTRICTED TO THE VEHICLE + CORE CONTEXT GENES
----------------------------------------------------------------------
The first version of this script scored Jaccard similarity over every IS and
cargo call in the +/-20 kb window. It found almost nothing: even loci known to
share the same vehicle at the same offset (the six VEHICLE-INTACT loci) scored
0.43-0.57, because the variable upstream region and the downstream tra tail
are -- by this project's own finding, drawn explicitly in the schematic panel
of colour_blocks.py -- "a different history in every locus". Including them
in a similarity score meant to find SHARED history actively works against it:
one locus-private IS anywhere in 40 kb costs real Jaccard when the feature set
is this small.

So the vocabulary here is the vehicle (--ends, default ISEsa1,ISEsa2) plus the
conserved core context genes (--core, default cytb,czc,relB,relE,hok_gef) --
the same five names draw_schematic() in colour_blocks.py already singles out
as "the conserved core". A fragment does not count as present: a 64 bp or
275 bp stub is evidence an event happened, not evidence the locus currently
carries a functional copy, and conflating the two would blur exactly the
VEHICLE-INTACT/VEHICLE-DECAYED distinction this grouping exists to sharpen.

GROUPING: UNION-FIND AT A SIMILARITY THRESHOLD
------------------------------------------------
Two loci are in the same group iff Jaccard(core signature A, core signature
B) >= --min-similarity. Groups are the connected components of that graph.
This can chain (A-B and B-C both clear the threshold without A-C clearing
it); --out flags every group where the weakest internal pair is below the
threshold so that is visible rather than hidden inside one colour. At the
default 0.90 on 14 loci it does not happen -- every multi-member group here
is pairwise 1.00 internally -- but it is checked, not assumed.

COLOUR IS CONTENT-ADDRESSED ON THE SHARED SIGNATURE, NOT ON MEMBERSHIP
-------------------------------------------------------------------------
Hashing the member list would reassign colour the moment a locus is added.
Hashing the sorted set of shared elements means a group keeps its colour
across runs as long as its architecture signature repeats -- the same
convention colour_blocks.py already uses for IS-family and gene-family
colour.
"""

import argparse, csv, glob, hashlib, itertools, sys
from collections import defaultdict

SPARE = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02',
         '#a6761d', '#666666', '#1f78b4', '#b2182b', '#542788', '#2d9f78',
         '#b15928', '#8c564b', '#17becf', '#bf5b17']


def group_colour(elements):
    if not elements:
        return '#cccccc'
    h = int(hashlib.sha1('|'.join(sorted(elements)).encode()).hexdigest(), 16)
    return SPARE[h % len(SPARE)]


def core_signature(rows, vocab):
    """Which vocabulary elements this locus carries as a complete or partial
    copy -- fragments excluded, see module docstring."""
    return {r['name'] for r in rows
            if r['name'] in vocab and r.get('cov_class') != 'fragment'}


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--blocks', required=True, help='directory of *.blocks.tsv')
    p.add_argument('--out', required=True)
    p.add_argument('--out-pairs', default='',
                   help='also write every pairwise Jaccard score here')
    p.add_argument('--ends', default='ISEsa1,ISEsa2',
                   help="vehicle elements -- must match the run's ENDS")
    p.add_argument('--core', default='cytb,czc,relB,relE,hok_gef',
                   help='conserved context genes, matching the schematic '
                        'panel in colour_blocks.py')
    p.add_argument('--min-similarity', type=float, default=0.90)
    a = p.parse_args()

    files = sorted(glob.glob(f'{a.blocks}/*.blocks.tsv'))
    if not files:
        sys.exit(f'no *.blocks.tsv under {a.blocks}')
    rows = {}
    for f in files:
        lid = f.split('/')[-1][:-len('.blocks.tsv')]
        rows[lid] = list(csv.DictReader(open(f), delimiter='\t'))

    vocab = {e.strip() for e in (a.ends + ',' + a.core).split(',') if e.strip()}
    loci = sorted(rows)
    sig = {l: core_signature(rows[l], vocab) for l in loci}

    pairs = {}
    for x, y in itertools.combinations(loci, 2):
        pairs[(x, y)] = jaccard(sig[x], sig[y])

    uf = UnionFind(loci)
    for (x, y), s in pairs.items():
        if s >= a.min_similarity:
            uf.union(x, y)

    members = defaultdict(list)
    for l in loci:
        members[uf.find(l)].append(l)

    groups = {}
    for root, ms in members.items():
        ms = sorted(ms)
        shared = (set.intersection(*(sig[m] for m in ms)) if len(ms) > 1
                  else sig[ms[0]])
        weakest = (min(pairs[tuple(sorted((x, y)))]
                      for x, y in itertools.combinations(ms, 2))
                   if len(ms) > 1 else 1.0)
        gid = 'G' + hashlib.sha1('|'.join(sorted(shared)).encode()).hexdigest()[:6]
        groups[root] = {'gid': gid, 'members': ms, 'shared': shared,
                        'weakest_pair': weakest, 'colour': group_colour(shared)}

    for g in groups.values():
        if len(g['members']) == 1:
            g['colour'] = '#cccccc'
            g['gid'] = 'SINGLETON'

    with open(a.out, 'w') as fh:
        fh.write('locus_id\tgroup_id\tgroup_colour\tgroup_size\t'
                  'group_members\tshared_core\tweakest_internal_pair\t'
                  'chained\n')
        for g in groups.values():
            chained = 'yes' if g['weakest_pair'] < a.min_similarity else 'no'
            for l in g['members']:
                fh.write('\t'.join([
                    l, g['gid'], g['colour'], str(len(g['members'])),
                    ','.join(g['members']), ','.join(sorted(g['shared'])) or '-',
                    f"{g['weakest_pair']:.2f}", chained]) + '\n')

    if a.out_pairs:
        with open(a.out_pairs, 'w') as fh:
            fh.write('locus_a\tlocus_b\tcore_jaccard\n')
            for (x, y), s in sorted(pairs.items(), key=lambda kv: -kv[1]):
                fh.write(f'{x}\t{y}\t{s:.3f}\n')

    multi = sorted((g for g in groups.values() if len(g['members']) > 1),
                   key=lambda g: -len(g['members']))
    singles = [g for g in groups.values() if len(g['members']) == 1]
    print(f'vocabulary: {", ".join(sorted(vocab))}')
    print(f'{len(loci)} loci -> {len(groups)} groups '
          f'(core Jaccard >= {a.min_similarity})\n')
    for g in multi:
        flag = ('  ** CHAINED -- not every internal pair clears the '
                'threshold **' if g['weakest_pair'] < a.min_similarity else '')
        print(f"  {g['gid']}  {len(g['members'])} loci, "
              f"weakest pair={g['weakest_pair']:.2f}{flag}")
        for l in g['members']:
            print(f"    {l.split('__')[0]}")
        print(f"    shared: {', '.join(sorted(g['shared'])) or '(none)'}\n")
    if singles:
        print(f"  singletons ({len(singles)}): "
              + ', '.join(s['members'][0].split('__')[0] for s in singles))
        for s in singles:
            print(f"    {s['members'][0].split('__')[0]}: "
                  f"{', '.join(sorted(s['shared'])) or '(none)'}")
    print(f'\n-> {a.out}')
    if a.out_pairs:
        print(f'-> {a.out_pairs}')


if __name__ == '__main__':
    main()
