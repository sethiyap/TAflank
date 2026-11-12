#!/usr/bin/env python3
"""
Shared sequence-homology block drawing/colouring library.

This used to be flanker_blocks.py -- a standalone script that built blocks
from pairwise dc-megablast and ALSO drew the figure. The blastn-based method
has been retired (pangraph is the sole block-finding engine now -- see
pangraph_blocks.py's own docstring for why), so the parts that were specific
to that method (load_hsps/chain_reach/build_blocks/UnionFind and the
argparse-driven main()) are gone. What's left is the part every remaining
caller (pangraph_blocks.py, pangraph_annotated.py, pangraph_context_summary.py)
actually imports: the colour/ordering primitives and the figure itself, which
don't care which engine produced the (locus -> [(start, end, member-set)])
blocks it's handed.

colour_for() gives the SAME colour to the SAME set of loci wherever it
appears (content-addressed via sha1), so "these N loci still share this much
flanking DNA" is visually legible without a legend. cluster_order() orders
figure rows by shared-sequence distance (hierarchical clustering + optimal
leaf ordering) so similar loci sit next to each other.
"""
import csv, hashlib, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_metadata import short_organism  # noqa: E402

CORE_COLOUR = '#f2a541'     # relB+relE itself -- fixed, not hashed
PRIVATE_COLOUR = '#fbe0ea'  # component of size 1: "nothing shared out here"
PRIVATE_EDGE = '#e7b9c9'


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def colour_for(members):
    if len(members) == 1:
        return PRIVATE_COLOUR, PRIVATE_EDGE
    h = hashlib.sha1(','.join(sorted(members)).encode()).hexdigest()
    hue = (int(h[:8], 16) % 1000) / 1000
    sat = 0.45 + (int(h[8:10], 16) % 30) / 100
    light = 0.52 + (int(h[10:12], 16) % 20) / 100
    import colorsys
    r, g, b = colorsys.hls_to_rgb(hue, light, sat)
    rgb = '#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255))
    return rgb, '#333333'


def cluster_order(loci, mutual, flank):
    """Leaf order from hierarchical clustering on pairwise sequence distance,
    so loci with the most shared flanking DNA sit next to each other in the
    figure instead of in accession order.

    distance(A,B) = 1 - (mutual_up + mutual_down) / (2*flank): 0 for a pair
    identical the whole window, towards 1 for a pair that shares nothing.
    Average-linkage + optimal_leaf_ordering, which (unlike a plain
    dendrogram) also minimises the distance between ADJACENT leaves
    specifically -- the thing that actually matters for "similar rows sit
    next to each other" rather than just "similar rows end up under the same
    branch somewhere".

    Falls back to the old span-based ordering if scipy is not installed,
    rather than failing the whole figure over a cosmetic improvement.
    """
    if len(loci) < 3:
        return sorted(loci)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.cluster.hierarchy import optimal_leaf_ordering
    except ImportError:
        print('WARNING: scipy not available -- falling back to span-based '
              'ordering (pip/conda install scipy for proper clustering)',
              file=sys.stderr)
        def span_score(l):
            return -sum(mutual.get((min(l, o), max(l, o)), {}).get(s, 0)
                       for o in loci if o != l for s in ('up', 'down'))
        return sorted(loci, key=lambda l: (span_score(l), l))

    n = len(loci)
    condensed = []
    for i in range(n):
        for j in range(i + 1, n):
            l1, l2 = loci[i], loci[j]
            m = mutual.get((min(l1, l2), max(l1, l2)), {'up': 0, 'down': 0})
            shared = m['up'] + m['down']
            condensed.append(1.0 - min(1.0, shared / (2 * flank)))
    Z = linkage(condensed, method='average')
    Z = optimal_leaf_ordering(Z, condensed)
    return [loci[i] for i in leaves_list(Z)]


def write_context_groups(path, loci, blocks_up, blocks_down, org_of):
    """Operates only on the blocks_up/blocks_down shape pangraph_blocks.py
    produces: {locus: [(start_d, end_d, frozenset(members))]}."""
    with open(path, 'w') as fh:
        fh.write('locus_id\torganism_short\tfurthest_shared_up\t'
                 'furthest_shared_down\tup_group_at_edge\t'
                 'down_group_at_edge\tfull_window_group\n')
        for l in loci:
            up_last_block = next((b for b in reversed(blocks_up[l])
                                  if len(b[2]) > 1), None)
            dn_last_block = next((b for b in reversed(blocks_down[l])
                                  if len(b[2]) > 1), None)
            up_edge = blocks_up[l][-1][2] if blocks_up[l] else frozenset([l])
            dn_edge = blocks_down[l][-1][2] if blocks_down[l] else frozenset([l])
            full = up_edge & dn_edge
            fh.write(f'{l}\t{org_of(l)}\t'
                     f'{up_last_block[1] if up_last_block else 0}\t'
                     f'{dn_last_block[1] if dn_last_block else 0}\t'
                     f'{",".join(sorted(up_edge))}\t'
                     f'{",".join(sorted(dn_edge))}\t'
                     f'{",".join(sorted(full))}\n')


def draw_figure(a, loci, order, blocks_up, blocks_down, ctx, groups, org_of,
                method_label='pangraph v1.4.0, N-way alignment graph',
                filename_stem='pangraph_blocks',
                context_groups_filename='pangraph_context_groups.tsv',
                row_scale=1.0):
    """Everything --no-figure skips. matplotlib is imported HERE, not at
    module level, so a cluster node with blast but no matplotlib can still
    run the caller's script for its TSV tables.

    `row_scale` halves (or otherwise scales) row height/pitch and the label
    font together -- see pangraph_context_summary.py's identical parameter
    for why font size must scale alongside the row rather than staying
    fixed (a fixed font overruns a half-height row)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # ---------------- figure ----------------
    PITCH, H = 1.0 * row_scale, 0.94 * row_scale
    label_fontsize = 7.5 * (0.5 + 0.5 * row_scale)
    fig_h = PITCH * len(loci) + 2.6
    fig, ax = plt.subplots(figsize=(14, fig_h))

    label_rows = []   # (locus_id, y, label_obj) -- swatch placed after
                       # axis limits are final, see the note below the loop
    for i, l in enumerate(order):
        y = (len(order) - 1 - i) * PITCH
        for side, bl in (('up', blocks_up[l]), ('down', blocks_down[l])):
            for d0, d1, mem in bl:
                x0 = -d1 if side == 'up' else d0
                x1 = -d0 if side == 'up' else d1
                if x1 <= x0:
                    continue
                fc, ec = colour_for(mem)
                ax.add_patch(Rectangle((x0, y - H / 2), x1 - x0, H,
                                       facecolor=fc, edgecolor='white',
                                       lw=0.4, zorder=2))
        ax.add_patch(Rectangle((a.core_start, y - H / 2),
                               a.core_end - a.core_start, H,
                               facecolor=CORE_COLOUR, edgecolor='white',
                               lw=0.5, zorder=3))
        acc = l.split('__')[0]
        org = org_of(l)
        vclass = ctx.get(l, {}).get('vehicle_class', '') if ctx else ''
        label_obj = ax.text(-a.flank * 1.03, y, f'{acc}  {org}', ha='right',
                            va='center', fontsize=label_fontsize,
                            style='italic' if org else 'normal')
        label_rows.append((l, y, label_obj))
        if vclass:
            ax.text(a.flank * 1.03, y, vclass, ha='left', va='center',
                    fontsize=label_fontsize - 1, color='#777777')

    ax.annotate('relB + relE\n(the anchor)', (a.core_end / 2, (len(order) - 1) * PITCH + 0.75),
               ha='center', va='bottom', fontsize=8, fontweight='bold',
               color='#b5620a',
               arrowprops=dict(arrowstyle='-', lw=0.8, color='#b5620a'),
               xytext=(a.core_end / 2, (len(order) - 1) * PITCH + 1.3))
    ax.axvline(0, color='#999999', lw=0.6, ls=':', zorder=1)

    ax.set_xlim(-a.flank * 1.18, a.flank * 1.12)
    ax.set_ylim(-0.8, (len(order) - 1) * PITCH + 1.9)

    # Group swatches, placed AFTER xlim/ylim are final and by MEASURING the
    # label text's own rendered width -- a gap sized as a fraction of the
    # x-range collapses under --xlim since the font size (and so the
    # label's rendered width) does not shrink with it.
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
            ax.add_patch(Rectangle((x_left, y - 0.17 * row_scale),
                                   x_right - x_left, 0.34 * row_scale,
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
    n_isolated = sum(
        1 for l in loci
        if all(len(m) == 1 for _, _, m in blocks_up[l]) and
           all(len(m) == 1 for _, _, m in blocks_down[l]))
    # Short, two-line title -- the long version (method description, full
    # legend-in-prose, pointer to the tsv) moved to this module's own
    # docstring and OUTPUTS.md; the figure itself only needs the number.
    ax.set_title(
        f'Sequence-homology blocks around relB+relE, {len(loci)} loci '
        f'({method_label})\n'
        f'{len(loci) - n_isolated}/{len(loci)} share flanking sequence with '
        '>=1 other locus; same colour = mutually alignable, pale pink = '
        'private',
        fontsize=10, pad=12)

    os.makedirs(a.outdir, exist_ok=True)
    fig.savefig(os.path.join(a.outdir, f'{filename_stem}.png'), dpi=200,
               bbox_inches='tight')
    fig.savefig(os.path.join(a.outdir, f'{filename_stem}.svg'),
               bbox_inches='tight')
    plt.close(fig)
