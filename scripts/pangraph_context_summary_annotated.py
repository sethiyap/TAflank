#!/usr/bin/env python3
"""
One row per UNIQUE context (pangraph_context_summary.py's model), but with
colour_blocks.py's gene/IS annotation arrows drawn on top of each row's
sequence background (pangraph_annotated.py's model) -- so the de-duplicated,
easy-to-scan row count of the summary figure and the "what actually sits
there" detail of the per-locus annotated figure don't have to be two
separate things to look at.

  pangraph_context_summary_annotated.py --graph results/pangraph/graph.json \
                                 --blocks run_all_v2/results/blocks \
                                 --groups results/pangraph_groups.tsv \
                                 --pairwise results/pangraph_pairwise.tsv \
                                 --context results/context.tsv \
                                 --junctions results/junctions.tsv \
                                 --summary results/summary.tsv \
                                 --outdir figures

Each row draws its CONTEXT's medoid locus (pangraph_context_summary.py's
medoid()) -- the member whose worst-case shared_frac to the rest of its
group is largest -- the same representative that script already picks, so
the two figures agree on which locus stands in for the group.
"""
import argparse, glob, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_metadata import short_organism  # noqa: E402
from colour_blocks import gene_arrow, simple_colour, SIMPLE, \
    load_junctions, dr_ticks  # noqa: E402
from seq_blocks import cluster_order, colour_for, CORE_COLOUR, \
    PRIVATE_COLOUR, read_tsv, savefig_capped  # noqa: E402
from pangraph_blocks import load_graph, locus_walk, core_range, \
    build_blocks_from_walk  # noqa: E402
from pangraph_context_summary import medoid  # noqa: E402

H_BG = 0.94
H_CDS, H_IS, H_CARGO = 0.10, 0.34, 0.34


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--graph', required=True)
    p.add_argument('--blocks', required=True,
                   help='dir of *.blocks.tsv (GENE_BLOCKS output)')
    p.add_argument('--groups', required=True,
                   help='results/pangraph_groups.tsv (groups.py)')
    p.add_argument('--pairwise', required=True,
                   help='results/pangraph_pairwise.tsv (pangraph_blocks.py)')
    p.add_argument('--context', default='')
    p.add_argument('--junctions', default='',
                   help='results/junctions.tsv (relbe_junctions.py)')
    p.add_argument('--summary', default='',
                   help='results/summary.tsv (collate.py) -- Inc-type badge')
    p.add_argument('--outdir', default='figures')
    p.add_argument('--flank', type=int, default=20000,
                   help='x-axis scale/ticks only, see pangraph_blocks.py')
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

    groups_rows = read_tsv(a.groups)
    if not groups_rows:
        sys.exit(f'no rows in {a.groups}')
    by_members, group_meta = {}, {}
    for r in groups_rows:
        key = frozenset(r['group_members'].split(','))
        by_members[key] = sorted(key)
        group_meta[key] = r

    pairwise_rows = read_tsv(a.pairwise)
    shared_frac = {(r['locus_a'], r['locus_b']): float(r['shared_frac'])
                  for r in pairwise_rows}
    mutual = {}
    for r in pairwise_rows:
        la, lb = r['locus_a'], r['locus_b']
        mutual[(min(la, lb), max(la, lb))] = {
            'up': float(r['shared_up_bp']), 'down': float(r['shared_down_bp'])}

    ctx = {}
    if a.context and os.path.exists(a.context):
        ctx = {r['locus_id']: r for r in read_tsv(a.context)}
    junctions = load_junctions(a.junctions)
    summary = {}
    if a.summary and os.path.exists(a.summary):
        summary = {r['locus_id']: r for r in read_tsv(a.summary)}

    def org_of(l):
        return short_organism(ctx.get(l, {}).get('organism', '') if ctx else '')

    keys = sorted(by_members, key=lambda k: (-len(by_members[k]), by_members[k]))
    contexts = []
    for i, key in enumerate(keys, 1):
        members = by_members[key]
        rep = medoid(members, shared_frac)
        contexts.append({
            'number': i, 'members': members, 'rep': rep,
            'colour': group_meta[key]['group_colour'],
        })

    feats_by_locus = {}
    for f in glob.glob(os.path.join(a.blocks, '*.blocks.tsv')):
        lid = os.path.basename(f)[:-len('.blocks.tsv')]
        if lid in walks:
            feats_by_locus[lid] = read_tsv(f)

    # Context NUMBER stays size-ranked (its stable identity); figure ROW
    # ORDER is cluster_order() on each context's medoid, using the same
    # shared-bp data that decides the blocks' own colours -- see
    # pangraph_context_summary.py's identical split and its longer comment
    # on why this is sequence similarity, not a vehicle_class/IS axis.
    rep_to_ctx = {c['rep']: c for c in contexts}
    draw_order = cluster_order(list(rep_to_ctx), mutual, a.flank)
    contexts_drawn = [rep_to_ctx[rep] for rep in draw_order]

    os.makedirs(a.outdir, exist_ok=True)

    PITCH = 1.0
    fig_h = PITCH * len(contexts) + 2.6
    fig, ax = plt.subplots(figsize=(14, fig_h))

    any_dr_drawn = False
    label_rows = []   # (y, text object, swatch colour)
    HEAD = a.flank / 90
    for i, c in enumerate(contexts_drawn):
        y = (len(contexts) - 1 - i) * PITCH
        rep = c['rep']

        # ---- background: pangraph blocks, representative locus only ----
        for side, bl in (('up', blocks_up[rep]), ('down', blocks_down[rep])):
            for d0, d1, mem in bl:
                x0 = -d1 if side == 'up' else d0
                x1 = -d0 if side == 'up' else d1
                if x1 <= x0:
                    continue
                fc, _ = colour_for(mem)
                ax.add_patch(Rectangle((x0, y - H_BG / 2), x1 - x0, H_BG,
                                       facecolor=fc, edgecolor='white',
                                       lw=0.4, zorder=1, alpha=0.9))
        ax.add_patch(Rectangle((core_lo, y - H_BG / 2), core_hi - core_lo,
                               H_BG, facecolor=CORE_COLOUR, edgecolor='white',
                               lw=0.5, zorder=2, alpha=0.9))

        # ---- foreground: gene/IS annotation arrows, representative locus --
        for r in feats_by_locus.get(rep, []):
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
            col = simple_colour(ft, nm, strand)
            frag = r.get('cov_class') == 'fragment'
            ax.add_patch(Polygon(
                gene_arrow(o, w, y, h, strand, HEAD), closed=True,
                facecolor=col, edgecolor='black', lw=0.7,
                hatch='///' if frag else None, zorder=4))
            dr = junctions.get((rep, nm, r['start'], r['end']))
            if dr:
                dr_ticks(ax, o, w, y, h, dr['same_dr_colour'])
                any_dr_drawn = True

        acc = rep.split('__')[0]
        n = len(c['members'])
        MAX_SHOWN = 8
        accs = [m.split('__')[0] for m in c['members']]
        shown = [f'{x}*' if n > 1 and x == acc else x for x in accs[:MAX_SHOWN]]
        mem_str = ', '.join(shown)
        if len(accs) > MAX_SHOWN:
            mem_str += f', +{len(accs) - MAX_SHOWN} more'
        rep_type = ctx.get(rep, {}).get('replicon_type', '') if ctx else ''
        inc = summary.get(rep, {}).get('plasmidfinder_replicons', '')
        rep_bits = [b for b in (rep_type, (inc if inc not in ('', '-') else ''))
                   if b]
        rep_txt = ('\n' + '  '.join(rep_bits)) if rep_bits else ''
        # Right-anchored just left of the plotted data -- text grows
        # LEFTWARD so it can never collide with the arrows to its right no
        # matter how long the member/Inc-type text gets. The swatch is
        # placed afterwards by MEASURING this text's rendered width (see
        # below the main loop), not a fixed data-unit offset guess -- a
        # fixed guess is what made an earlier version of this figure (and
        # pangraph_context_summary.py's) draw the swatch and label on top
        # of each other, and here additionally on top of the arrows.
        label_obj = ax.text(-a.flank * 1.02, y,
                f"Context {c['number']:02d}  (n={n})  {org_of(rep)}\n"
                f"{mem_str}{rep_txt}",
                ha='right', va='center', fontsize=6.6, linespacing=1.3)
        label_rows.append((y, label_obj, c['colour']))

    ax.annotate('relB + relE\n(the anchor)',
               ((core_lo + core_hi) / 2, (len(contexts) - 1) * PITCH + 0.75),
               ha='center', va='bottom', fontsize=8, fontweight='bold',
               color='#b5620a',
               arrowprops=dict(arrowstyle='-', lw=0.8, color='#b5620a'),
               xytext=((core_lo + core_hi) / 2, (len(contexts) - 1) * PITCH + 1.3))
    ax.axvline(0, color='#999999', lw=0.6, ls=':', zorder=0)
    ax.set_xlim(-a.flank * 1.46, a.flank * 1.12)
    ax.set_ylim(-0.8, (len(contexts) - 1) * PITCH + 1.9)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for y, label_obj, colour in label_rows:
        bbox = label_obj.get_window_extent(renderer=renderer)
        y_px = ax.transData.transform((0, y))[1]
        inv = ax.transData.inverted()
        (x_right, _) = inv.transform((bbox.x0 - 6, y_px))   # 6 px gap
        x_left = x_right - a.flank * 0.045                  # fixed swatch width
        ax.add_patch(Rectangle((x_left, y - 0.26), x_right - x_left, 0.52,
                               facecolor=colour, edgecolor='black',
                               lw=0.5, clip_on=False))

    step = 5000
    ticks = list(range(-a.flank, a.flank + 1, step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{v/1000:+.0f} kb' if v else '0' for v in ticks],
                       fontsize=8)
    ax.set_yticks([])
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel('distance from relB start (bp)', fontsize=9)
    dr_note = ('\n▽▽ at an element\'s edges = exact target-site duplication '
              'shared with another locus, colour = which exact repeat'
              if any_dr_drawn else '')
    ax.set_title(
        f'{len(contexts)} unique contexts summarising {len(loci)} loci -- '
        'annotated (pangraph background + colour_blocks.py arrows)\n'
        'one row = one context, drawn from its medoid (most-typical) '
        'member; left-edge swatch = context-group colour (pangraph_groups.tsv)\n'
        'arrows = named IS/cargo/CDS; strand is the way each arrow points, '
        'not its colour\n'
        'row order = shared sequence (same signal as the block colours) -- '
        'numbers out of sequence sit next to the one they share DNA with'
        + dr_note,
        fontsize=9.5, pad=16)

    legend_items = [
        ('sequence: shared', dict(facecolor='#6baed6')),
        ('sequence: private', dict(facecolor=PRIVATE_COLOUR)),
        ('sequence: core (relB+relE)', dict(facecolor=CORE_COLOUR)),
        ('IS (arrow = strand)', dict(facecolor=SIMPLE['is'])),
        ('relB', dict(facecolor=SIMPLE['relB'])),
        ('relE', dict(facecolor=SIMPLE['relE'])),
        ('cytb', dict(facecolor=SIMPLE['cytb'])),
        ('czc', dict(facecolor=SIMPLE['czc'])),
        ('other cargo', dict(facecolor=SIMPLE['cargo'])),
    ]
    for i, (label, spec) in enumerate(legend_items):
        yy = -1.4 - (i // 5) * 0.4
        xx = -a.flank * 1.44 + (i % 5) * a.flank * 0.5
        ax.add_patch(Rectangle((xx, yy - 0.12), a.flank * 0.03, 0.24,
                               edgecolor='black', lw=0.4, clip_on=False,
                               **spec))
        ax.text(xx + a.flank * 0.045, yy, label, fontsize=6.5, va='center',
               clip_on=False)

    savefig_capped(fig, os.path.join(a.outdir, 'pangraph_context_summary_annotated.png'),
                  dpi=200, bbox_inches='tight')
    savefig_capped(fig, os.path.join(a.outdir, 'pangraph_context_summary_annotated.svg'),
                  bbox_inches='tight')
    plt.close(fig)
    print(f'{len(contexts)} contexts ({len(loci)} loci) -> '
          f'{a.outdir}/pangraph_context_summary_annotated.png  .svg')


if __name__ == '__main__':
    main()
