#!/usr/bin/env python3
"""
Gene/IS annotation overlay (colour_blocks.py's arrows), drawn on top of
pangraph_blocks.py's sequence-homology background.

  pangraph_annotated.py --graph results/pangraph/graph.json \
                         --blocks run_all_v2/results/blocks \
                         --context run_all_v2/results/context.tsv \
                         --groups run_all_v2/results/pangraph_groups.tsv \
                         --outdir figures

The annotation layer (what gene/IS sits where) comes from GENE_BLOCKS'
output -- blastn-vs-ISfinder/cargo-search, plus Bakta's gene calls when
`params.run_bakta` is on (gene_blocks.py's annotation source order is Bakta
GFF, then parent GenBank, then none) -- independent of which engine finds
the SEQUENCE-homology blocks underneath it.
"""
import argparse, glob, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_metadata import short_organism  # noqa: E402
from colour_blocks import gene_arrow, simple_colour, SIMPLE  # noqa: E402
from seq_blocks import cluster_order, colour_for, CORE_COLOUR, \
    PRIVATE_COLOUR, read_tsv  # noqa: E402
from pangraph_blocks import load_graph, locus_walk, core_range, \
    build_blocks_from_walk, mutual_from_blocks  # noqa: E402

H_BG = 0.94
H_CDS, H_IS, H_CARGO = 0.10, 0.34, 0.34


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--graph', required=True)
    p.add_argument('--blocks', required=True,
                   help='dir of *.blocks.tsv (GENE_BLOCKS output)')
    p.add_argument('--context', default='')
    p.add_argument('--groups', default='')
    p.add_argument('--outdir', default='figures')
    p.add_argument('--flank', type=int, default=20000,
                   help='x-axis scale/ticks only, see pangraph_blocks.py')
    a = p.parse_args()

    graph = load_graph(a.graph)
    walks, flank_up = locus_walk(graph)
    loci = sorted(walks)
    if not loci:
        sys.exit(f'no paths in {a.graph}')

    from collections import defaultdict
    block_members = defaultdict(set)
    for lid, entries in walks.items():
        for bid, _, _, _ in entries:
            block_members[bid].add(lid)
    block_members = {k: frozenset(v) for k, v in block_members.items()}

    core_lo, core_hi = core_range(walks)
    blocks_up, blocks_down = build_blocks_from_walk(
        walks, block_members, core_lo, core_hi)

    # Derived from blocks_up/blocks_down (already trimmed to the core
    # boundary), not from raw walks -- see pangraph_blocks.py's
    # mutual_from_blocks() docstring for why that distinction matters.
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

    feats_by_locus = {}
    for f in glob.glob(os.path.join(a.blocks, '*.blocks.tsv')):
        lid = os.path.basename(f)[:-len('.blocks.tsv')]
        if lid in walks:
            feats_by_locus[lid] = read_tsv(f)

    os.makedirs(a.outdir, exist_ok=True)

    PITCH = 1.0
    fig_h = 1.0 * PITCH * len(loci) + 2.6
    fig, ax = plt.subplots(figsize=(14, fig_h))

    label_rows = []
    for i, l in enumerate(order):
        y = (len(order) - 1 - i) * PITCH

        # ---- background: pangraph blocks ----
        for side, bl in (('up', blocks_up[l]), ('down', blocks_down[l])):
            for d0, d1, mem in bl:
                x0 = -d1 if side == 'up' else d0
                x1 = -d0 if side == 'up' else d1
                if x1 <= x0:
                    continue
                fc, ec = colour_for(mem)
                ax.add_patch(Rectangle((x0, y - H_BG / 2), x1 - x0, H_BG,
                                       facecolor=fc, edgecolor='white',
                                       lw=0.4, zorder=1, alpha=0.9))
        ax.add_patch(Rectangle((core_lo, y - H_BG / 2), core_hi - core_lo,
                               H_BG, facecolor=CORE_COLOUR, edgecolor='white',
                               lw=0.5, zorder=2, alpha=0.9))

        # ---- foreground: gene/IS annotation arrows (colour_blocks.py) ----
        HEAD = a.flank / 90
        for r in feats_by_locus.get(l, []):
            ft, nm = r['feature_type'], r['name']
            o = int(r['offset_from_bait'])
            w = int(r['end']) - int(r['start']) + 1
            vo, ve = max(o, -a.flank), min(o + w, a.flank)
            if ve <= vo:
                continue
            o, w = vo, ve - vo
            strand = r.get('strand', '+')
            if ft == 'CDS':
                ax.add_patch(Polygon(
                    gene_arrow(o, w, y, H_CDS, strand, HEAD), closed=True,
                    facecolor=SIMPLE['cds'], edgecolor='#999999',
                    lw=0.25, zorder=3, alpha=0.85))
                continue
            h = H_IS if ft == 'IS' else H_CARGO
            c = simple_colour(ft, nm, strand)
            frag = r.get('cov_class') == 'fragment'
            ax.add_patch(Polygon(
                gene_arrow(o, w, y, h, strand, HEAD), closed=True,
                facecolor=c, edgecolor='black', lw=0.7,
                hatch='///' if frag else None, zorder=4))

        acc = l.split('__')[0]
        org = org_of(l)
        label_obj = ax.text(-a.flank * 1.03, y, f'{acc}  {org}', ha='right',
                            va='center', fontsize=7.5,
                            style='italic' if org else 'normal')
        label_rows.append((l, y, label_obj))
        vclass = ctx.get(l, {}).get('vehicle_class', '') if ctx else ''
        if vclass:
            ax.text(a.flank * 1.03, y, vclass, ha='left', va='center',
                    fontsize=6.5, color='#777777')

    ax.annotate('relB + relE\n(the anchor)',
               ((core_lo + core_hi) / 2, (len(order) - 1) * PITCH + 0.75),
               ha='center', va='bottom', fontsize=8, fontweight='bold',
               color='#b5620a',
               arrowprops=dict(arrowstyle='-', lw=0.8, color='#b5620a'),
               xytext=((core_lo + core_hi) / 2, (len(order) - 1) * PITCH + 1.3))
    ax.axvline(0, color='#999999', lw=0.6, ls=':', zorder=0)

    ax.set_xlim(-a.flank * 1.18, a.flank * 1.12)
    ax.set_ylim(-0.8, (len(order) - 1) * PITCH + 1.9)

    if groups:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for lid, y, label_obj in label_rows:
            grp = groups.get(lid)
            if not grp or grp['group_colour'] == '#cccccc':
                continue
            bbox = label_obj.get_window_extent(renderer=renderer)
            y_px = ax.transData.transform((0, y))[1]
            inv = ax.transData.inverted()
            (x_right, _) = inv.transform((bbox.x0 - 5, y_px))
            (x_left, _) = inv.transform((bbox.x0 - 15, y_px))
            ax.add_patch(Rectangle((x_left, y - 0.17), x_right - x_left, 0.34,
                                   facecolor=grp['group_colour'],
                                   edgecolor='none', clip_on=False))

    step = 5000
    ticks = list(range(-a.flank, a.flank + 1, step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{v/1000:+.0f} kb' if v else '0' for v in ticks],
                       fontsize=8)
    ax.set_yticks([])
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel('distance from relB start (bp)', fontsize=9)
    ax.set_title(
        f'relBE flanking architecture: annotation on top of pangraph '
        f'sequence identity, {len(loci)} loci\n'
        'background colour = pangraph block (Noll, Molari & Neher 2023); '
        'pale pink = private, nothing shared beyond this point\n'
        'arrows = named IS/cargo/CDS (colour_blocks.py, "simple" style) -- '
        'the SAME arrow colour on a DIFFERENT background is the same named '
        'element on different DNA',
        fontsize=9.5, pad=14)

    legend_items = [
        ('sequence: shared', dict(facecolor='#6baed6')),
        ('sequence: private', dict(facecolor=PRIVATE_COLOUR)),
        ('sequence: core (relB+relE)', dict(facecolor=CORE_COLOUR)),
        ('IS, forward', dict(facecolor=SIMPLE['is_fwd'])),
        ('IS, reverse', dict(facecolor=SIMPLE['is_rev'])),
        ('relB', dict(facecolor=SIMPLE['relB'])),
        ('relE', dict(facecolor=SIMPLE['relE'])),
        ('cytb', dict(facecolor=SIMPLE['cytb'])),
        ('czc', dict(facecolor=SIMPLE['czc'])),
        ('other cargo', dict(facecolor=SIMPLE['cargo'])),
    ]
    for i, (label, spec) in enumerate(legend_items):
        yy = -1.4 - (i // 5) * 0.4
        xx = -a.flank * 1.1 + (i % 5) * a.flank * 0.5
        ax.add_patch(Rectangle((xx, yy - 0.12), a.flank * 0.03, 0.24,
                               edgecolor='black', lw=0.4, clip_on=False,
                               **spec))
        ax.text(xx + a.flank * 0.045, yy, label, fontsize=6.5, va='center',
               clip_on=False)

    fig.savefig(os.path.join(a.outdir, 'pangraph_annotated.png'), dpi=200,
               bbox_inches='tight')
    fig.savefig(os.path.join(a.outdir, 'pangraph_annotated.svg'),
               bbox_inches='tight')
    plt.close(fig)
    print(f'{len(loci)} loci -> {a.outdir}/pangraph_annotated.png  .svg')


if __name__ == '__main__':
    main()
