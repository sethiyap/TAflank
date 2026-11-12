#!/usr/bin/env python3
"""
Stage 6 input -- group loci by shared SEQUENCE, the direct counterpart to
structural_groups.py's gene-identity Jaccard grouping.

Generic: takes any pairwise table with locus_a/locus_b/shared_frac columns
-- currently that's pangraph_pairwise.tsv (pangraph_blocks.py), used by the
PANGRAPH_GROUPS process in main.nf.

  groups.py --pairwise results/pangraph_pairwise.tsv \
            --out results/pangraph_groups.tsv

Mirrors structural_groups.py's method exactly -- union-find at a similarity
threshold, group colour content-addressed on the group's own identity rather
than its membership list -- but the similarity is `shared_frac` (fraction of
the 2x flank window two loci are still mutually alignable, per pangraph's
alignment graph) instead of Jaccard over IS/cargo gene names.

This table is a CLASSIFICATION output, read by pairwise_evidence.py-style
reconciliation and drawn as a swatch in pangraph_blocks.py's figure. It does
NOT drive that figure's row order -- ordering stays the continuous-distance
clustering in seq_blocks.py's cluster_order().
"""
import argparse, csv, hashlib, itertools, sys
from collections import defaultdict

# Same palette shape as structural_groups.py's SPARE, deliberately DIFFERENT
# literal colours -- a sequence-based group and a structural (gene-name)
# group are different claims (same DNA vs. same named genes) and must never
# be visually confusable even when, as often happens, their membership
# coincides.
SPARE = ['#2166ac', '#b2182b', '#4daf4a', '#984ea3', '#ff7f00', '#a65628',
         '#f781bf', '#999999', '#1b9e77', '#e6ab02', '#66c2a5', '#fc8d62',
         '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']


def group_colour(key):
    if not key:
        return '#cccccc'
    h = int(hashlib.sha1(key.encode()).hexdigest(), 16)
    return SPARE[h % len(SPARE)]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


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
    p.add_argument('--pairwise', required=True,
                   help='results/pangraph_pairwise.tsv (pangraph_blocks.py)')
    p.add_argument('--out', required=True)
    p.add_argument('--min-similarity', type=float, default=0.85,
                   help='shared_frac >= this = same group (union-find edge)')
    a = p.parse_args()

    rows = read_tsv(a.pairwise)
    if not rows:
        sys.exit(f'no rows in {a.pairwise}')
    loci = sorted({r['locus_a'] for r in rows} | {r['locus_b'] for r in rows})
    pairs = {(r['locus_a'], r['locus_b']): float(r['shared_frac'])
             for r in rows}
    org = {}
    for r in rows:
        org.setdefault(r['locus_a'], r.get('organism_short_a', ''))
        org.setdefault(r['locus_b'], r.get('organism_short_b', ''))

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
        pair_scores = ([pairs[tuple(sorted((x, y)))]
                        for x, y in itertools.combinations(ms, 2)]
                       if len(ms) > 1 else [1.0])
        weakest = min(pair_scores)
        # content-addressed on the GROUP's own identity (its sorted member
        # list), same spirit as structural_groups.py hashing the shared
        # signature -- membership is this group's only "identity", since
        # unlike the gene-based grouping there is no small named vocabulary
        # to hash instead
        gid = 'G' + hashlib.sha1('|'.join(ms).encode()).hexdigest()[:6]
        groups[root] = {'gid': gid, 'members': ms, 'weakest_pair': weakest,
                        'colour': group_colour('|'.join(ms))}

    for g in groups.values():
        if len(g['members']) == 1:
            g['colour'] = '#cccccc'
            g['gid'] = 'SINGLETON'

    with open(a.out, 'w') as fh:
        fh.write('locus_id\torganism_short\tgroup_id\tgroup_colour\t'
                 'group_size\tgroup_members\tweakest_internal_pair\tchained\n')
        for g in groups.values():
            chained = 'yes' if g['weakest_pair'] < a.min_similarity else 'no'
            # weakest_internal_pair is the LEAST shared_frac between any two
            # members of this group -- same meaning as
            # structural_groups.py's column of the same name, just computed
            # from sequence shared_frac instead of gene-name Jaccard. A group
            # never implies more similarity than its weakest link supports.
            for l in g['members']:
                fh.write('\t'.join([
                    l, org.get(l, ''), g['gid'], g['colour'],
                    str(len(g['members'])), ','.join(g['members']),
                    f"{g['weakest_pair']:.3f}", chained]) + '\n')

    multi = sorted((g for g in groups.values() if len(g['members']) > 1),
                   key=lambda g: -len(g['members']))
    singles = [g for g in groups.values() if len(g['members']) == 1]
    print(f'{len(loci)} loci -> {len(groups)} groups '
          f'(shared sequence fraction >= {a.min_similarity})\n')
    for g in multi:
        flag = ('  ** CHAINED -- not every internal pair clears the '
                'threshold **' if g['weakest_pair'] < a.min_similarity else '')
        print(f"  {g['gid']}  {len(g['members'])} loci, "
              f"weakest pair={g['weakest_pair']:.3f}{flag}")
        for l in g['members']:
            print(f"    {l.split('__')[0]}")
    if singles:
        print(f"\n  singletons ({len(singles)}): "
              + ', '.join(s['members'][0].split('__')[0] for s in singles))
    print(f'\n-> {a.out}')


if __name__ == '__main__':
    main()
