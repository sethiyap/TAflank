#!/usr/bin/env python3
"""
Stage 7 -- the Acman Fig 4-style colour-block figure, aligned on the bait.

  colour_blocks.py --blocks results/blocks --context results/context.tsv \
                   --metadata results/record_metadata.tsv --outdir figures

Produces
  colour_blocks.png/.svg   one row per locus, x = offset from the bait
  is_matrix.png            loci x IS element presence/absence
  is_matrix.tsv            the same matrix as a table (Stage 6 input)

WHY THE CDS TRACK IS GREY AND NOT COLOURED
------------------------------------------
Acman's middle panel colours every gene by its pangenome family. We cannot do
that yet and must not pretend to: 153 of the 278 annotation-named CDS in this
set are "hypothetical protein", so `colour_key` collapses them onto ONE value.
Colouring by it would draw a single huge block of one colour across every locus
and invite the reading that those regions are homologous, which is exactly the
claim the figure would be made to support. Until MMseqs2 clustering exists
(Stage 3), CDS are drawn as uncoloured context and only the two tracks that
carry trustworthy names are coloured:

  IS     -- ISfinder element names, from is_search.py
  cargo  -- the curated gene set, from cargo_search.py

COLOUR STABILITY
----------------
Cluster indices renumber whenever the dataset grows, which would recolour the
figure between runs and break comparison with anything already published. So
the biologically load-bearing names are hand-pinned and everything else is
content-addressed by a hash of its own name -- an element keeps its colour when
new loci are added.
"""

import argparse, csv, glob, hashlib, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.lines import Line2D

# Hand-pinned: these carry the argument, so their colours must never drift.
PINNED = {
    'ISEsa1': '#d62728',   # the candidate ancestral vehicle
    'ISEsa2': '#ff7f0e',   # its partner at the -3898 boundary
    'relB':   '#1f77b4',   # the bait
    'relE':   '#17becf',
    'czc':    '#2ca02c',   # cobalt-zinc-cadmium, the metal-resistance core
    'cytb':   '#98df8a',
    'hok_gef': '#9467bd',
    'tra':    '#8c564b',
    'intI':   '#e377c2',
    'ars':    '#7f7f7f',
    'sil':    '#bcbd22',
    'qnr':    '#c5b0d5',
}
SPARE = ['#aec7e8', '#ffbb78', '#c7c7c7', '#dbdb8d', '#f7b6d2',
         '#c49c94', '#9edae5', '#ff9896', '#c5b0d5', '#8dd3c7',
         '#fdb462', '#b3de69', '#fccde5', '#bc80bd', '#ccebc5']


def colour(name):
    if name in PINNED:
        return PINNED[name]
    h = int(hashlib.sha1(name.encode()).hexdigest(), 16)
    return SPARE[h % len(SPARE)]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--blocks', required=True, help='directory of *.blocks.tsv')
    p.add_argument('--context', required=True)
    p.add_argument('--metadata', default='')
    p.add_argument('--outdir', default='figures')
    p.add_argument('--flank', type=int, default=20000)
    a = p.parse_args()

    files = sorted(glob.glob(os.path.join(a.blocks, '*.blocks.tsv')))
    if not files:
        sys.exit(f'no *.blocks.tsv under {a.blocks}')
    rows = []
    for f in files:
        rows += read_tsv(f)

    ctx = {r['locus_id']: r for r in read_tsv(a.context)}
    meta = {}
    if a.metadata and os.path.exists(a.metadata):
        for r in read_tsv(a.metadata):
            meta[r['accession']] = r

    by_locus = {}
    for r in rows:
        by_locus.setdefault(r['locus_id'], []).append(r)

    # Order loci by vehicle class, then by the ISEsa1 offset -- so loci that
    # share an architecture sit next to each other and a shared offset reads as
    # a vertical line rather than needing to be hunted for.
    CLASS_ORDER = {'COMPLETE': 0, 'DECAYED': 1, 'HALF-CORE': 2, 'NO CORE': 3}

    def esa1_offset(lid):
        offs = [int(r['offset_from_bait']) for r in by_locus[lid]
                if r['feature_type'] == 'IS' and r['name'] == 'ISEsa1']
        return min(offs) if offs else 1 << 30

    loci = sorted(by_locus, key=lambda l: (
        CLASS_ORDER.get(ctx.get(l, {}).get('vehicle_class', ''), 9),
        esa1_offset(l), l))

    # ---------------- panel A: aligned colour blocks ----------------
    PITCH = 1.75          # vertical distance between loci
    fig, ax = plt.subplots(figsize=(15, 0.95 * PITCH * len(loci) + 2.6))
    H_CDS, H_IS, H_CARGO = 0.17, 0.40, 0.30
    seen = {}

    for i, lid in enumerate(loci):
        y = (len(loci) - 1 - i) * PITCH
        label_y = []      # x-ranges already labelled on this row, for stacking
        feats = by_locus[lid]
        offs = [int(r['offset_from_bait']) for r in feats]
        lo, hi = min(offs + [-a.flank]), max(offs + [a.flank])
        ax.plot([lo, hi], [y, y], color='#d9d9d9', lw=1.0, zorder=0,
                solid_capstyle='butt')

        for r in feats:
            o = int(r['offset_from_bait'])
            w = int(r['end']) - int(r['start']) + 1
            ft, nm = r['feature_type'], r['name']
            clipped = r.get('clipped_at_window_edge', 'no') not in ('no', '-', '')

            if ft == 'CDS':
                # context only: no family colours exist yet (see docstring)
                ax.add_patch(Rectangle((o, y - H_CDS / 2), w, H_CDS,
                                       facecolor='#f0f0f0', edgecolor='#bdbdbd',
                                       lw=0.3, zorder=1))
                continue

            h = H_IS if ft == 'IS' else H_CARGO
            yy = y + 0.34 if ft == 'IS' else y - 0.36
            c = colour(nm)
            frag = r.get('cov_class') == 'fragment'
            ax.add_patch(Rectangle(
                (o, yy - h / 2), w, h, facecolor=c,
                edgecolor='black' if not clipped else '#555555',
                lw=0.9 if not clipped else 0.9,
                linestyle='-' if not clipped else (0, (1.4, 1.0)),
                hatch='///' if frag else None, alpha=0.95, zorder=3))
            seen.setdefault(nm, (c, ft))
            # name the vehicle and the small-but-decisive fragments
            if nm.startswith('ISEsa') or (frag and ft == 'IS'):
                # Short elements sit shoulder to shoulder and their labels
                # overprint into an unreadable smear. Push each one up a level
                # until it clears every label already placed on this row.
                cx, half = o + w / 2, max(w, 1300) / 2
                lvl = 0
                while any(abs(cx - px) < (half + ph) and lv == lvl
                          for px, ph, lv in label_y):
                    lvl += 1
                label_y.append((cx, half, lvl))
                ax.annotate(f'{nm}\n{w} bp' if frag else nm,
                            (cx, yy + h / 2 + 0.04 + 0.30 * lvl),
                            ha='center', va='bottom', fontsize=5.6,
                            color='#222222', zorder=5, linespacing=0.95,
                            arrowprops=dict(arrowstyle='-', lw=0.4,
                                            color='#999999') if lvl else None)

        if not any(r['feature_type'] == 'CDS' for r in feats):
            ax.text(-a.flank * 1.02, y - 0.62, 'no CDS in the GenBank record',
                    fontsize=5.8, style='italic', color='#aa0000',
                    ha='left', va='center')
        cls = ctx.get(lid, {}).get('vehicle_class', '?')
        acc = lid.split('__')[0]
        org = meta.get(acc, {}).get('organism', '')
        org = ' '.join(org.split()[:2]) if org else ''
        ax.text(-a.flank * 1.045, y, f'{acc}\n{org}', ha='right', va='center',
                fontsize=7.2, linespacing=1.25, style='italic' if org else 'normal')
        ax.text(a.flank * 1.02, y, cls, ha='left', va='center', fontsize=7,
                color={'COMPLETE': '#1a7f1a', 'DECAYED': '#b8860b',
                       'HALF-CORE': '#d2691e'}.get(cls, '#999999'),
                fontweight='bold')

    ax.axvline(0, color='#1f77b4', lw=1.3, ls='--', zorder=4)
    ax.text(0, (len(loci) - 1) * PITCH + 0.95, ' relB start (bait = 0)', fontsize=7.5,
            color='#1f77b4', ha='left', va='bottom')
    # the recurring boundary, drawn so it can be checked rather than asserted
    ax.axvline(-3898, color='#ff7f0e', lw=0.9, ls=':', zorder=2)
    ax.text(-3898, -1.35, '-3898', fontsize=6.5, color='#ff7f0e', ha='center')

    ax.set_xlim(-a.flank * 1.30, a.flank * 1.22)
    ax.set_ylim(-1.55, (len(loci) - 1) * PITCH + 1.30)
    ax.set_yticks([])
    ax.set_xlabel('position relative to relB start (bp)', fontsize=9)
    ax.set_xticks(range(-a.flank, a.flank + 1, 5000))
    ax.set_xticklabels([f'{v//1000:+d} kb' if v else '0'
                        for v in range(-a.flank, a.flank + 1, 5000)], fontsize=7.5)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)

    def sq(c, lab, ec='black'):
        return Line2D([0], [0], marker='s', ls='', markersize=7,
                      markerfacecolor=c, markeredgecolor=ec,
                      markeredgewidth=0.5, label=lab)

    is_names = sorted(n for n in seen if seen[n][1] == 'IS')
    cg_names = sorted(n for n in seen if seen[n][1] != 'IS')
    handles = ([sq('none', 'IS elements:', 'none')]
               + [sq(seen[n][0], n) for n in is_names]
               + [sq('none', '', 'none')] * ((-len(is_names) - 1) % 7)
               + [sq('none', 'cargo genes:', 'none')]
               + [sq(seen[n][0], n) for n in cg_names])
    handles += [
        Line2D([0], [0], marker='s', ls='', markersize=7, markerfacecolor='#f0f0f0',
               markeredgecolor='#bdbdbd', label='CDS (no family colour yet)'),
        Line2D([0], [0], marker='s', ls='', markersize=7, markerfacecolor='white',
               markeredgecolor='black', label='hatched = fragment (<50% cov)'),
    ]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.11),
              ncol=7, fontsize=6.6, frameon=False, handletextpad=0.35,
              columnspacing=1.0)
    ax.set_title('relBE flanking architecture, aligned on the bait\n'
                 'IS track above each line, curated cargo below; CDS grey '
                 '(families not yet called)', fontsize=10, pad=12)

    os.makedirs(a.outdir, exist_ok=True)
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(a.outdir, f'colour_blocks.{ext}'),
                    dpi=200, bbox_inches='tight')
    plt.close(fig)

    # ---------------- panel B: IS presence/absence heatmap ----------------
    is_rows = [r for r in rows if r['feature_type'] == 'IS']
    elems = {}
    for r in is_rows:
        elems.setdefault(r['name'], []).append(int(r['offset_from_bait']))
    # order by how many loci carry it, then by median offset: the recurring
    # elements collect on the left, the one-off graffiti trails off right
    n_loci_with = {e: len({r['locus_id'] for r in is_rows if r['name'] == e})
                   for e in elems}
    order_e = sorted(elems, key=lambda e: (-n_loci_with[e],
                                           sorted(elems[e])[len(elems[e]) // 2]))

    cell = {}
    for r in is_rows:
        k = (r['locus_id'], r['name'])
        prev = cell.get(k)
        if prev is None or r['cov_class'] == 'complete':
            cell[k] = r['cov_class']

    figb, axb = plt.subplots(figsize=(0.42 * len(order_e) + 4.2,
                                      0.46 * len(loci) + 2.0))
    CMAP = {'complete': '#08519c', 'partial': '#6baed6', 'fragment': '#c6dbef'}
    for i, lid in enumerate(loci):
        y = len(loci) - 1 - i
        for j, e in enumerate(order_e):
            cv = cell.get((lid, e))
            axb.add_patch(Rectangle((j, y), 1, 1,
                                    facecolor=CMAP.get(cv, '#f7f7f7'),
                                    edgecolor='white', lw=1.1))
    axb.set_xlim(0, len(order_e)); axb.set_ylim(0, len(loci))
    axb.set_xticks([j + 0.5 for j in range(len(order_e))])
    axb.set_xticklabels(order_e, rotation=90, fontsize=7)
    axb.set_yticks([len(loci) - 1 - i + 0.5 for i in range(len(loci))])
    axb.set_yticklabels([f"{l.split('__')[0]}  "
                         f"{ctx.get(l, {}).get('vehicle_class', '?')}"
                         for l in loci], fontsize=7)
    axb.set_title('IS element presence by locus\n'
                  f'ordered by number of loci carrying it (n={len(loci)})',
                  fontsize=9.5)
    for s in ('top', 'right', 'left', 'bottom'):
        axb.spines[s].set_visible(False)
    axb.tick_params(length=0)
    axb.legend(handles=[Line2D([0], [0], marker='s', ls='', markersize=8,
                               markerfacecolor=v, markeredgecolor='white',
                               label=k) for k, v in CMAP.items()],
               loc='upper center', bbox_to_anchor=(0.5, -0.32), ncol=3,
               fontsize=7, frameon=False)
    figb.tight_layout()
    figb.savefig(os.path.join(a.outdir, 'is_matrix.png'), dpi=200,
                 bbox_inches='tight')
    plt.close(figb)

    with open(os.path.join(a.outdir, 'is_matrix.tsv'), 'w') as fh:
        fh.write('locus_id\tvehicle_class\t' + '\t'.join(order_e) + '\n')
        for lid in loci:
            fh.write(lid + '\t' + ctx.get(lid, {}).get('vehicle_class', '?'))
            for e in order_e:
                fh.write('\t' + (cell.get((lid, e)) or '0'))
            fh.write('\n')

    print(f'{len(loci)} loci, {len(is_rows)} IS calls, '
          f'{len(order_e)} distinct elements')
    print(f'  recurring (>=2 loci): '
          f'{", ".join(e for e in order_e if n_loci_with[e] >= 2)}')
    print(f'  single-locus        : '
          f'{sum(1 for e in order_e if n_loci_with[e] == 1)} elements')
    print(f'-> {a.outdir}/colour_blocks.png  .svg')
    print(f'-> {a.outdir}/is_matrix.png  is_matrix.tsv')


if __name__ == '__main__':
    main()
