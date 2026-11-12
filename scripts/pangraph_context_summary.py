#!/usr/bin/env python3
"""
One row per UNIQUE context (a pangraph_groups.tsv group), numbered, plus a
two-sheet Excel workbook listing every locus in each context -- built from
pangraph's blocks.

  pangraph_context_summary.py --graph results/pangraph/graph.json \
                               --groups results/pangraph_groups.tsv \
                               --pairwise results/pangraph_pairwise.tsv \
                               --context results/context.tsv \
                               --outdir figures --out-xlsx results/pangraph_context_summary.xlsx

Medoid-representative logic: the member whose worst-case shared_frac to the
rest of its group is largest is drawn as that context's one representative
row.
"""
import argparse, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_metadata import short_organism  # noqa: E402
from seq_blocks import colour_for, CORE_COLOUR, PRIVATE_COLOUR, read_tsv  # noqa: E402
from pangraph_blocks import load_graph, locus_walk, core_range, \
    build_blocks_from_walk  # noqa: E402


def medoid(members, pairs):
    """Same rule as context_summary.py's: the member whose worst (minimum)
    shared_frac to any other member in the group is the LARGEST."""
    if len(members) == 1:
        return members[0]
    best, best_score = None, -1.0
    for m in sorted(members):
        worst = min(pairs.get(tuple(sorted((m, o))), 0.0)
                   for o in members if o != m)
        if worst > best_score:
            best, best_score = m, worst
    return best


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--graph', required=True)
    p.add_argument('--groups', required=True,
                   help='results/pangraph_groups.tsv (groups.py '
                        '--pairwise pangraph_pairwise.tsv)')
    p.add_argument('--pairwise', required=True,
                   help='results/pangraph_pairwise.tsv (pangraph_blocks.py)')
    p.add_argument('--context', default='')
    p.add_argument('--outdir', default='figures')
    p.add_argument('--out-xlsx', required=True)
    p.add_argument('--flank', type=int, default=20000,
                   help='x-axis scale/ticks only, see pangraph_blocks.py')
    p.add_argument('--row-scale', type=float, default=0.5,
                   help='row HEIGHT multiplier, see context_summary.py\'s '
                        'own --row-scale for the full explanation')
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
    # Key contexts by their ACTUAL member set, not group_id -- see
    # context_summary.py's identical comment: every singleton shares the
    # literal group_id 'SINGLETON'.
    by_members = {}
    group_meta = {}
    for r in groups_rows:
        key = frozenset(r['group_members'].split(','))
        by_members[key] = sorted(key)
        group_meta[key] = r

    pairwise_rows = read_tsv(a.pairwise)
    shared_frac = {(r['locus_a'], r['locus_b']): float(r['shared_frac'])
                  for r in pairwise_rows}

    ctx = {}
    if a.context and os.path.exists(a.context):
        ctx = {r['locus_id']: r for r in read_tsv(a.context)}

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
            'weakest': group_meta[key].get('weakest_internal_pair', '1.000'),
        })

    # ---------------- Excel workbook ----------------
    import openpyxl
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = openpyxl.Workbook()
    ws1 = wb.active
    ws1.title = 'Contexts'
    hdr1 = ['context_number', 'n_loci', 'representative_locus',
           'representative_accession', 'representative_organism',
           'weakest_internal_shared_frac', 'members']
    ws1.append(hdr1)
    for c in ws1[1]:
        c.font = Font(bold=True)
    for c in contexts:
        ws1.append([
            c['number'], len(c['members']), c['rep'],
            c['rep'].split('__')[0], org_of(c['rep']),
            c['weakest'], ', '.join(m.split('__')[0] for m in c['members']),
        ])
        ws1.cell(row=ws1.max_row, column=1).fill = PatternFill(
            'solid', fgColor=c['colour'].lstrip('#'))
    for i, w in enumerate((14, 8, 34, 14, 22, 16, 50), 1):
        ws1.column_dimensions[get_column_letter(i)].width = w

    ws2 = wb.create_sheet('Loci')
    ctx_hdr = ['accession', 'organism', 'genus', 'replicon_type',
              'vehicle_class', 'core_class', 'cargo_families', 'n_is']
    ws2.append(['context_number', 'is_representative', 'locus_id'] + ctx_hdr)
    for c in ws2[1]:
        c.font = Font(bold=True)
    loci_to_ctx = {m: c for c in contexts for m in c['members']}
    for l in loci:
        c = loci_to_ctx.get(l)
        row = ctx.get(l, {})
        ws2.append([
            c['number'] if c else '', 'yes' if c and c['rep'] == l else 'no', l,
        ] + [row.get(k, '') for k in ctx_hdr])
        if c:
            ws2.cell(row=ws2.max_row, column=1).fill = PatternFill(
                'solid', fgColor=c['colour'].lstrip('#'))
    for i, w in enumerate((14, 15, 34, 12, 34, 14, 14, 14, 24, 6), 1):
        ws2.column_dimensions[get_column_letter(i)].width = w
    ws2.freeze_panes = 'A2'
    ws1.freeze_panes = 'A2'

    os.makedirs(os.path.dirname(a.out_xlsx) or '.', exist_ok=True)
    wb.save(a.out_xlsx)
    print(f'{len(contexts)} contexts ({len(loci)} loci) -> {a.out_xlsx}')

    # ---------------- figure: one row per context ----------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # See context_summary.py's identical lines for why label_fontsize must
    # scale alongside PITCH/H rather than staying fixed.
    PITCH, H = 1.0 * a.row_scale, 0.94 * a.row_scale
    label_fontsize = 7.2 * (0.5 + 0.5 * a.row_scale)
    fig_h = PITCH * len(contexts) + 2.2
    fig, ax = plt.subplots(figsize=(14, fig_h))

    for i, c in enumerate(contexts):
        y = (len(contexts) - 1 - i) * PITCH
        rep = c['rep']
        for side, bl in (('up', blocks_up[rep]), ('down', blocks_down[rep])):
            for d0, d1, mem in bl:
                x0 = -d1 if side == 'up' else d0
                x1 = -d0 if side == 'up' else d1
                if x1 <= x0:
                    continue
                fc, _ = colour_for(mem)
                ax.add_patch(Rectangle((x0, y - H / 2), x1 - x0, H,
                                       facecolor=fc, edgecolor='white',
                                       lw=0.4, zorder=2))
        ax.add_patch(Rectangle((core_lo, y - H / 2), core_hi - core_lo, H,
                               facecolor=CORE_COLOUR, edgecolor='white',
                               lw=0.5, zorder=3))
        ax.add_patch(Rectangle((-a.flank * 1.30, y - 0.26 * a.row_scale),
                               a.flank * 0.045, 0.52 * a.row_scale,
                               facecolor=c['colour'], edgecolor='black',
                               lw=0.5, clip_on=False))
        acc = rep.split('__')[0]
        n = len(c['members'])
        MAX_SHOWN = 10
        accs = [m.split('__')[0] for m in c['members']]
        shown = [f'{x}*' if n > 1 and x == acc else x for x in accs[:MAX_SHOWN]]
        mem_str = ', '.join(shown)
        if len(accs) > MAX_SHOWN:
            mem_str += f', +{len(accs) - MAX_SHOWN} more'
        ax.text(-a.flank * 1.27, y,
                f"Context {c['number']:02d}  (n={n})\n{mem_str}",
                ha='left', va='center', fontsize=label_fontsize, linespacing=1.3)

    ax.annotate('relB + relE\n(the anchor)',
               ((core_lo + core_hi) / 2, (len(contexts) - 1) * PITCH + 0.75),
               ha='center', va='bottom', fontsize=8, fontweight='bold',
               color='#b5620a',
               arrowprops=dict(arrowstyle='-', lw=0.8, color='#b5620a'),
               xytext=((core_lo + core_hi) / 2, (len(contexts) - 1) * PITCH + 1.3))
    ax.axvline(0, color='#999999', lw=0.6, ls=':', zorder=1)
    ax.set_xlim(-a.flank * 1.32, a.flank * 1.12)
    ax.set_ylim(-0.8, (len(contexts) - 1) * PITCH + 1.9)
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
        f'{len(contexts)} unique contexts summarising {len(loci)} loci '
        '(pangraph)\n'
        'one row = one context, drawn from its medoid (most-typical) '
        f'member; see {os.path.basename(a.out_xlsx)} for every locus '
        'in each context\n'
        'block colour = pangraph sequence identity (same colour = same set '
        'of loci sharing that exact block); pale pink = private to this '
        'context; left-edge swatch = context-group colour',
        fontsize=10, pad=16)

    legend_items = [
        ('sequence: shared (colour = which loci share it)',
         dict(facecolor='#6baed6')),
        ('sequence: private (unique to this context)',
         dict(facecolor=PRIVATE_COLOUR)),
        ('sequence: core (relB + relE)', dict(facecolor=CORE_COLOUR)),
        ('left-edge swatch = context group (pangraph_groups.tsv), '
         'not a sequence block', dict(facecolor='#999999')),
    ]
    for i, (label, spec) in enumerate(legend_items):
        yy = -1.4 - i * 0.4
        xx = -a.flank * 1.30
        ax.add_patch(Rectangle((xx, yy - 0.12), a.flank * 0.04, 0.24,
                               edgecolor='black', lw=0.4, clip_on=False,
                               **spec))
        ax.text(xx + a.flank * 0.06, yy, label, fontsize=6.5, va='center',
               clip_on=False)

    os.makedirs(a.outdir, exist_ok=True)
    fig.savefig(os.path.join(a.outdir, 'pangraph_context_summary.png'),
               dpi=200, bbox_inches='tight')
    fig.savefig(os.path.join(a.outdir, 'pangraph_context_summary.svg'),
               bbox_inches='tight')
    plt.close(fig)
    print(f'-> {a.outdir}/pangraph_context_summary.png  .svg')


if __name__ == '__main__':
    main()
