#!/usr/bin/env python3
"""
Builds the sequence-homology blocks around relB+relE (which loci share which
flanking DNA, drawn as colour blocks) from a pangraph alignment graph.

  pangraph_blocks.py --graph results/pangraph/graph.json \
                     --context results/context.tsv --groups results/pangraph_groups.tsv \
                     --outdir figures --out-tables results

pangraph (Noll, Molari & Neher 2023, Microbial Genomics 9:001034;
github.com/neherlab/pangraph) builds a true N-way alignment graph: every
input genome is an explicit ordered walk over blocks, each block with its
own strand, so a rearrangement (inversion, internal reordering) is
represented directly rather than approximated by chaining pairwise hits
outward from an anchor and hoping nothing flipped along the way.

This script does NOT re-run the paralog/completeness filter (bait_search.py)
-- it reads --graph, already built (by build_pangraph_graph.py) from the
filtered, extracted loci windows.

REUSE, NOT REIMPLEMENTATION
------------------------------------------------------------------------
`draw_figure`, `cluster_order`, `colour_for`, `CORE_COLOUR`,
`PRIVATE_COLOUR` are imported from seq_blocks.py unchanged. Only the thing
that actually differs between a block-finding engine and the next -- how a
locus's offset axis gets cut into (start, end, member-set) blocks in the
first place -- is written fresh here; everything downstream of that (colour
assignment, row ordering, the figure itself, the table schemas) is shared
code in seq_blocks.py.
"""
import argparse, json, os, re, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_metadata import short_organism  # noqa: E402
from seq_blocks import (draw_figure, cluster_order, colour_for,
                        CORE_COLOUR, PRIVATE_COLOUR, read_tsv,
                        write_context_groups)  # noqa: E402

FLANK_UP_RE = re.compile(r'flank_up=(\d+)')


def load_graph(path):
    with open(path) as fh:
        return json.load(fh)


def locus_walk(graph):
    """-> {locus_id: [(block_id, offset_start, offset_end, strand), ...]},
    {locus_id: flank_up} -- offsets are relative to relB start (0), using
    EACH locus's own flank_up from EXTRACT_LOCI's own fasta description,
    not an assumed-uniform CLI value -- correct even for a locus whose
    window was clamped at a contig edge, where flank_up legitimately
    differs from the rest."""
    nodes = graph['nodes']
    walks, flank_up = {}, {}
    for p in graph['paths'].values():
        lid = p['name']
        m = FLANK_UP_RE.search(p.get('desc', ''))
        if not m:
            sys.exit(f'{lid}: no flank_up= in the fasta description -- '
                     f'was this window built by extract_loci.py?')
        fu = int(m.group(1))
        flank_up[lid] = fu
        entries = []
        for nid in p['nodes']:
            n = nodes[str(nid)]
            s, e = n['position']
            entries.append((n['block_id'], s - fu, e - fu, n['strand']))
        entries.sort(key=lambda t: t[1])
        walks[lid] = entries
    return walks, flank_up


def core_range(walks):
    """The core block's offset span per locus (whichever block covers
    offset 0, i.e. relB start), then the UNION across all loci -- so every
    row draws the core band at the same coordinates, and a locus whose own
    core-block alignment is a few bp narrower than another's doesn't make
    its neighbouring blocks misclassified as core/non-core."""
    spans = {}
    for lid, entries in walks.items():
        for bid, o0, o1, strand in entries:
            if o0 <= 0 < o1:
                spans[lid] = (o0, o1)
                break
        else:
            sys.exit(f'{lid}: no block covers offset 0 (relB start) -- '
                     f'check --graph was built from windows with relB at '
                     f'the same place in every fasta')
    return min(s for s, _ in spans.values()), max(e for _, e in spans.values())


def build_blocks_from_walk(walks, block_members, core_lo, core_hi):
    """-> blocks_up, blocks_down in the shape draw_figure() (imported from
    seq_blocks.py, unmodified) expects: {locus: [(d0, d1,
    frozenset(members))]}. draw_figure plots these as
    `x0 = -d1 if side=='up' else d0; x1 = -d0 if side=='up' else d1` -- NOT
    offset by core_lo/core_hi -- so for that to land each block at its true
    absolute position, d0/d1 must be -o1/-o0 on the up side (plain negation,
    no core_lo term) and o0/o1 directly on the down side (no core_hi term).
    An earlier version of this function subtracted core_lo/core_hi here,
    which is what the earlier blastn-chaining method's own block-builder
    effectively did too (its 'down' distances were computed relative to a
    fixed core edge and plotted unadjusted) -- invisible there only because
    its core span (527 bp, relB+relE alone) was small enough that the fixed
    CORE_COLOUR rectangle (drawn at a higher zorder) happened to paint over
    the resulting misalignment. Pangraph's real core span is wider
    (core_hi-core_lo, not a fixed 527), so that masking can't be relied on
    here -- hence the direct absolute conversion below instead.

    Boundaries are not inferred by raising a threshold over pairwise reach,
    either -- they are each locus's OWN real block boundaries from the
    graph, which is exactly pangraph's advantage: no approximation step.
    """
    up, down = defaultdict(list), defaultdict(list)
    for lid, entries in walks.items():
        for bid, o0, o1, strand in entries:
            if o1 <= core_lo:
                up[lid].append((-o1, -o0, block_members[bid]))
            elif o0 >= core_hi:
                down[lid].append((o0, o1, block_members[bid]))
            else:
                # Straddles the GLOBAL core window -- core_lo/core_hi is the
                # union across all loci of each one's own core-covering
                # node, so a locus whose own core ends a little earlier
                # than that union can have its very next (genuinely
                # downstream) block start a few bp inside [core_lo,
                # core_hi]. Dropping the whole node here -- as an earlier
                # version of this function did -- silently discards that
                # downstream content instead of just the sliver that
                # actually overlaps the core; for a 14/40000bp pair this
                # dropped nearly an ENTIRE 14kb shared downstream block,
                # giving a shared_frac of 0.35 against the true ~0.88.
                # Trim to the non-core remainder(s) instead -- the portion
                # still inside [core_lo, core_hi] is masked by the fixed
                # CORE_COLOUR rectangle regardless, so dropping only that
                # sliver changes nothing visible.
                if o0 < core_lo:
                    up[lid].append((-core_lo, -o0, block_members[bid]))
                if o1 > core_hi:
                    down[lid].append((core_hi, o1, block_members[bid]))
    for d in (up, down):
        for lid in d:
            d[lid].sort(key=lambda t: t[0])
            merged = []
            for d0, d1, mem in d[lid]:
                if merged and merged[-1][2] == mem and merged[-1][1] == d0:
                    merged[-1] = (merged[-1][0], d1, mem)
                else:
                    merged.append((d0, d1, mem))
            d[lid] = merged
    return up, down


def mutual_from_blocks(loci, blocks_up, blocks_down):
    """-> {(l1,l2): {'up': shared_bp, 'down': shared_bp}}, derived from the
    ALREADY-TRIMMED blocks_up/blocks_down (build_blocks_from_walk's output)
    rather than reclassifying raw walk entries against core_lo/core_hi a
    second time -- an earlier version did the latter and silently dropped
    any block straddling the (per-locus-variable) core boundary, which for
    one real pair understated shared_frac as 0.35 against a true ~0.88.
    Deriving it from blocks_up/down instead guarantees this number is
    consistent with what the figure actually draws, by construction."""
    mutual = {l1l2: {'up': 0, 'down': 0}
             for i, l1 in enumerate(loci) for l2 in loci[i + 1:]
             for l1l2 in [(l1, l2)]}
    for side, blocks in (('up', blocks_up), ('down', blocks_down)):
        for l1 in loci:
            for d0, d1, mem in blocks[l1]:
                w = d1 - d0
                for l2 in mem:
                    # l1's own blocks_up/down list AND l2's both contain
                    # this same shared block (block_members[bid] is
                    # symmetric), so counting every l2 in mem here would
                    # double-count each pair once from each side -- only
                    # add when l2 > l1, i.e. exactly once per unordered
                    # pair, using l1's own locally-measured width as that
                    # pair's shared width (loci's own copies of "the same"
                    # block differ by at most a few bp from indels)
                    if l2 <= l1:
                        continue
                    mutual[(l1, l2)][side] += w
    return mutual


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--graph', required=True,
                   help='graph.json from build_pangraph_graph.py')
    p.add_argument('--context', default='')
    p.add_argument('--groups', default='',
                   help='a *_groups.tsv in groups.py\'s schema -- swatch '
                        'only, does not change row order')
    p.add_argument('--outdir', default='figures')
    p.add_argument('--out-tables', default='results')
    p.add_argument('--flank', type=int, default=20000,
                   help='for the x-axis scale/ticks only -- real offsets '
                        'come from each locus\'s own flank_up metadata')
    p.add_argument('--row-scale', type=float, default=0.5,
                   help='row HEIGHT (and label font) multiplier -- 0.5 '
                        '(default) packs rows tight so many loci still fit '
                        'in one figure; 1.0 is the old, more spread-out '
                        'look')
    p.add_argument('--no-figure', action='store_true')
    a = p.parse_args()

    graph = load_graph(a.graph)
    walks, flank_up = locus_walk(graph)
    loci = sorted(walks)
    if not loci:
        sys.exit(f'no paths in {a.graph}')

    block_members = defaultdict(set)
    for lid, entries in walks.items():
        for bid, _, _, _ in entries:
            block_members[bid].add(lid)
    block_members = {k: frozenset(v) for k, v in block_members.items()}

    core_lo, core_hi = core_range(walks)
    blocks_up, blocks_down = build_blocks_from_walk(
        walks, block_members, core_lo, core_hi)

    # mutual[(l1,l2)] = {'up': shared_bp, 'down': shared_bp} -- feeds
    # cluster_order() (row ordering), derived from blocks_up/blocks_down
    # (already trimmed to the core boundary), not recomputed from raw
    # walks -- see mutual_from_blocks()'s docstring for why that matters.
    mutual = mutual_from_blocks(loci, blocks_up, blocks_down)

    ctx = {}
    if a.context and os.path.exists(a.context):
        ctx = {r['locus_id']: r for r in read_tsv(a.context)}
    groups = {}
    if a.groups and os.path.exists(a.groups):
        groups = {r['locus_id']: r for r in read_tsv(a.groups)}

    def org_of(l):
        return short_organism(ctx.get(l, {}).get('organism', '') if ctx else '')

    order = cluster_order(loci, mutual, a.flank)

    os.makedirs(a.out_tables, exist_ok=True)
    with open(os.path.join(a.out_tables, 'pangraph_pairwise.tsv'), 'w') as fh:
        # Same column schema groups.py (generic Jaccard-style grouping)
        # already expects, so PANGRAPH_GROUPS is just
        # `groups.py --pairwise pangraph_pairwise.tsv`, no dedicated script
        fh.write('locus_a\tlocus_b\torganism_short_a\torganism_short_b\t'
                 'shared_up_bp\tshared_down_bp\tshared_total_bp\tshared_frac\n')
        denom = 2 * max(flank_up.values())
        for (l1, l2), m in sorted(mutual.items()):
            tot = m['up'] + m['down']
            fh.write(f'{l1}\t{l2}\t{org_of(l1)}\t{org_of(l2)}\t'
                     f'{m["up"]}\t{m["down"]}\t{tot}\t{tot / denom:.4f}\n')

    with open(os.path.join(a.out_tables, 'pangraph_blocks.tsv'), 'w') as fh:
        fh.write('locus_id\torganism_short\tside\tstart_offset\tend_offset\t'
                 'group_size\tgroup_members\n')
        for l in loci:
            for side, bl in (('up', blocks_up[l]), ('down', blocks_down[l])):
                for d0, d1, mem in bl:
                    # d0/d1 are already true absolute offsets from relB
                    # start (see build_blocks_from_walk's docstring)
                    s0 = -d1 if side == 'up' else d0
                    s1 = -d0 if side == 'up' else d1
                    fh.write(f'{l}\t{org_of(l)}\t{side}\t{s0}\t{s1}\t'
                             f'{len(mem)}\t{",".join(sorted(mem))}\n')

    write_context_groups(os.path.join(a.out_tables, 'pangraph_context_groups.tsv'),
                        loci, blocks_up, blocks_down, org_of)

    if not a.no_figure:
        class NS:
            pass
        fa = NS()
        fa.flank, fa.core_start, fa.core_end = a.flank, core_lo, core_hi
        fa.outdir, fa.out_tables = a.outdir, a.out_tables
        draw_figure(fa, loci, order, blocks_up, blocks_down, ctx, groups,
                   org_of, method_label='pangraph v1.4.0, N-way alignment graph',
                   filename_stem='pangraph_blocks',
                   context_groups_filename='pangraph_context_groups.tsv',
                   row_scale=a.row_scale)

    n_private = sum(1 for l in loci
                    if all(len(m) == 1 for _, _, m in blocks_up[l]) and
                       all(len(m) == 1 for _, _, m in blocks_down[l]))
    print(f'{len(loci)} loci, {len(block_members)} distinct blocks, '
          f'core span {core_lo}..{core_hi} ({core_hi - core_lo} bp)')
    print(f'{len(loci) - n_private}/{len(loci)} loci share some flanking '
          f'sequence with at least one other locus')
    print(f'-> {a.out_tables}/pangraph_blocks.tsv')
    print(f'-> {a.out_tables}/pangraph_pairwise.tsv')
    print(f'-> {a.out_tables}/pangraph_context_groups.tsv')
    if not a.no_figure:
        print(f'-> {a.outdir}/pangraph_blocks.png  .svg')


if __name__ == '__main__':
    main()
