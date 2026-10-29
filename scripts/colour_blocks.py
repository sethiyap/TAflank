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

COLOUR: FAMILY IS THE HUE, ELEMENT IS THE SHADE
-----------------------------------------------
Twenty-four elements given twenty-four unrelated colours is a lookup table, not
a picture -- nothing can be read off it without crossing to the legend. IS names
already carry their family (ISfinder `ISEsa1_ISNCY_IS1202` -> family ISNCY), so
the family sets the hue and elements within it differ only in lightness. Twelve
hues instead of twenty-four, and related elements now look related.

Two things this makes visible that per-element colour hid:

  * ISEsa1 is **ISNCY** and ISEsa2 is **IS110** -- the two ends of the candidate
    vehicle are NOT the same family, so the "two IS of the same family flanking
    the cargo" composite-transposon rule does not describe it.
  * ISEsa1 and ISKpn21 **are** the same family and group (ISNCY/IS1202). They
    now share a hue, which is precisely the pair that makes a naive
    same-family rule fire on a nested element that flanks nothing.

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
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D

import colorsys

# ---------------------------------------------------------------------------
# Palettes lifted from the project's own decks, by reading the shape fills out
# of the .pptx XML rather than sampling a screenshot.
#
#   deck     relBE_schematic.pptx slide 11 -- the nine-loci figure this
#            pipeline reproduces, and the source CLAUDE.md names. Conserved
#            core pale, ISEsa1-like its own teal, variable IS warm.
#   lab2018  RelBE_from_Sept2018_lab_meeting.pptx slide 2. Highlights only the
#            TA pair: relB green, relE red, everything else white.
#
# In both, the hues that carry meaning are pinned exactly; families the deck
# never drew are given tints constrained to the deck's own warm range, so the
# upstream still reads as "variable IS" at a glance while staying separable
# on a closer look.
# ---------------------------------------------------------------------------
PALETTES = {
    'deck': {
        'cargo': {'relB': '#028090', 'relE': '#d1495b',
                  'cytb': '#d8e0e0', 'czc': '#d8e0e0',
                  'qnr': '#b57f45', 'hok_gef': '#7e8c8d', 'tra': '#4a5859',
                  'ars': '#8a8a84', 'sil': '#b6c2c2', 'intI': '#e09f3e'},
        'family': {'ISNCY': '#5ba4a8', 'IS110': '#b57f45'},
        'warm': (0.055, 0.115),      # hue range for every other IS family
        'cds_sat': 0.42,
    },
    'lab2018': {
        'cargo': {'relB': '#00b050', 'relE': '#ff0000',
                  'cytb': '#f2f2f2', 'czc': '#f2f2f2', 'qnr': '#f2f2f2',
                  'hok_gef': '#a5a5a5', 'tra': '#7f7f7f', 'ars': '#bfbfbf',
                  'sil': '#d9d9d9', 'intI': '#ed7d31'},
        'family': {'ISNCY': '#00b0f0', 'IS110': '#ffc000'},
        'warm': (0.055, 0.115),
        'cds_sat': 0.42,
    },
    'default': None,                 # the hand-built palette below
}
ACTIVE = {}

# Family hues. ISNCY and IS110 keep the red/orange of the elements that made
# them matter, so the figure stays comparable with earlier versions.
FAMILY = {
    'ISNCY':     '#d62728',   # ISEsa1, ISKpn21
    'IS110':     '#ff7f0e',   # ISEsa2
    'IS3':       '#2b7bba',
    'IS5':       '#8c6d31',
    'IS6':       '#6a3d9a',
    'IS1':       '#17becf',
    'IS4':       '#b15928',
    'IS21':      '#1b9e77',
    'IS91':      '#e7298a',
    'IS481':     '#66a61e',
    'IS1380':    '#7570b3',
    'IS200/IS605': '#a6761d',
}
# Darkest shade inside a family goes to the element that carries the argument.
ELEMENT_RANK = {'ISEsa1': 0, 'ISEsa2': 0, 'ISKpn21': 1}


def shade(base, i, n):
    """i-th of n lightness steps on `base`; i=0 is the base colour itself."""
    if n <= 1 or i == 0:
        return base
    r, g, b = (int(base[k:k + 2], 16) / 255 for k in (1, 3, 5))
    h, l, sat = colorsys.rgb_to_hls(r, g, b)
    l = min(0.90, l + (0.86 - l) * (i / n))
    sat = max(0.18, sat * (1 - 0.30 * i / n))
    return '#%02x%02x%02x' % tuple(
        int(round(c * 255)) for c in colorsys.hls_to_rgb(h, l, sat))


def gene_arrow(x, w, y, h, strand, head):
    """Block pointed in the direction of transcription.

    Vertical position in this figure is the TRACK (IS above the line, CDS on
    it, cargo below) and says nothing about strand -- which was read as strand
    the first time someone looked at it. So strand is drawn explicitly, and
    the only place it lives is the direction of the point.

    `head` is capped at half the block so a short feature stays a block with a
    point rather than collapsing into a bare triangle.
    """
    hd = min(head, w * 0.5)
    if strand == '-':
        return [(x, y), (x + hd, y + h / 2), (x + w, y + h / 2),
                (x + w, y - h / 2), (x + hd, y - h / 2)]
    return [(x, y - h / 2), (x + w - hd, y - h / 2), (x + w, y),
            (x + w - hd, y + h / 2), (x, y + h / 2)]


def cds_family_colour(fid):
    """Pastel hue from the content-addressed family id.

    Gene families are not named in a legend -- there are far too many, and the
    figure does not ask you to identify one. It asks whether the SAME colour
    appears at the same offset in another locus. So the only requirement is
    that a family's colour is stable and distinguishable from its neighbours,
    which a hash of the id gives for free. Pastel, so the CDS track stays
    visually behind the IS and cargo blocks that carry the argument.
    """
    hue = (int(fid[1:9], 16) % 3600) / 3600
    return '#%02x%02x%02x' % tuple(
        int(round(c * 255))
        for c in colorsys.hls_to_rgb(hue, 0.78, ACTIVE.get('cds_sat', 0.55)))


def family_colour(fam):
    pin = ACTIVE.get('family')
    if pin is not None:
        if fam in pin:
            return pin[fam]
        # constrained to the deck's warm range so every IS the deck did not
        # name still looks like the deck's "variable upstream", not confetti
        lo, hi = ACTIVE['warm']
        h = int(hashlib.sha1(fam.encode()).hexdigest(), 16)
        hue = lo + (h % 1000) / 1000 * (hi - lo)
        light = 0.38 + ((h >> 10) % 100) / 100 * 0.22
        return '#%02x%02x%02x' % tuple(
            int(round(c * 255))
            for c in colorsys.hls_to_rgb(hue, light, 0.55))
    if fam in FAMILY:
        return FAMILY[fam]
    h = int(hashlib.sha1(fam.encode()).hexdigest(), 16)
    return SPARE[h % len(SPARE)]


# Hand-pinned cargo colours: these carry the argument, so they must not drift.
PINNED = {
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
    cargo = ACTIVE.get('cargo')
    if cargo and name in cargo:
        return cargo[name]
    if name in PINNED:
        return PINNED[name]
    h = int(hashlib.sha1(name.encode()).hexdigest(), 16)
    return SPARE[h % len(SPARE)]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))



def legend_blocks(fam_seen, is_colour, seen):
    """Legend content as uniform (header, colour, [(label, swatch-spec)]) blocks.

    Everything -- families, cargo, the hatch/dash key -- is one kind of object,
    so the layout code below does not special-case any of them."""
    fams = sorted(fam_seen, key=lambda f: (-len(fam_seen[f]), f))
    blocks = [(f, family_colour(f),
               [(e, dict(facecolor=is_colour[e]))
                for e in sorted(fam_seen[f],
                                key=lambda x: (ELEMENT_RANK.get(x, 9), x))])
              for f in fams]
    cargo = sorted(n for n in seen if seen[n][1] != 'IS')
    for i in range(0, len(cargo), 6):
        blocks.append(('cargo genes' if i == 0 else '', '#333333',
                       [(n, dict(facecolor=seen[n][0])) for n in cargo[i:i + 6]]))
    blocks.append(('key', '#333333', [
        ('CDS: shared family', dict(facecolor='#b9d8ea', edgecolor='#9e9e9e')),
        ('CDS: this locus only', dict(facecolor='#f0f0f0', edgecolor='#bdbdbd')),
        ('fragment <50% cov', dict(facecolor='white', hatch='///')),
        ('runs off the window', dict(facecolor='white',
                                     linestyle=(0, (1.4, 1.0)))),
        ('point = strand', dict(facecolor='#dddddd')),
        ('row position = track', dict(facecolor='none', edgecolor='none')),
    ]))
    return blocks


def draw_legend(lg, blocks, ncol=6):
    """Pack the blocks into `ncol` columns, each row as tall as its tallest."""
    rows, row = [], []
    for b in blocks:
        row.append(b)
        if len(row) == ncol:
            rows.append(row); row = []
    if row:
        rows.append(row)
    heights = [max(len(b[2]) for b in r) + 1.45 for r in rows]
    total = sum(heights)
    lg.set_xlim(0, ncol); lg.set_ylim(total, 0); lg.axis('off')
    y0 = 0.0
    for r, h in zip(rows, heights):
        for k, (name, hcol, items) in enumerate(r):
            if name:
                lg.text(k + 0.03, y0 + 0.5, name, fontsize=7.2,
                        fontweight='bold', va='center', color=hcol)
            for j, (label, spec) in enumerate(items):
                yy = y0 + 1.45 + j
                lg.add_patch(Rectangle(
                    (k + 0.06, yy - 0.30), 0.115, 0.56,
                    edgecolor=spec.pop('edgecolor', 'black'), lw=0.55,
                    clip_on=False, **spec))
                lg.text(k + 0.24, yy, label, fontsize=6.5, va='center')
        y0 += h
    return total


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--blocks', required=True, help='directory of *.blocks.tsv')
    p.add_argument('--context', required=True)
    p.add_argument('--metadata', default='')
    p.add_argument('--outdir', default='figures')
    p.add_argument('--flank', type=int, default=20000)
    p.add_argument('--families', default='',
                   help='results/gene_families.tsv -- colours the CDS track')
    p.add_argument('--palette', choices=sorted(PALETTES), default='deck',
                   help="'deck' = relBE_schematic.pptx slide 11 (default), "
                        "'lab2018' = the Sept 2018 lab-meeting deck, "
                        "'default' = this script's own palette")
    p.add_argument('--min-family-loci', type=int, default=2,
                   help='colour a gene family only if it spans this many loci '
                        '(1 would give every singleton its own colour and turn '
                        'the track into confetti)')
    a = p.parse_args()

    ACTIVE.clear()
    if PALETTES[a.palette]:
        ACTIVE.update(PALETTES[a.palette])

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

    # element -> family, and a colour per element: family hue, element shade
    fam_of = {r['name']: (r['family'] or 'unclassified')
              for r in rows if r['feature_type'] == 'IS'}
    fam_members = {}
    for e, f in fam_of.items():
        fam_members.setdefault(f, []).append(e)
    is_colour = {}
    for f, members in fam_members.items():
        members.sort(key=lambda e: (ELEMENT_RANK.get(e, 9), e))
        base = family_colour(f)
        for i, e in enumerate(members):
            is_colour[e] = shade(base, i, len(members))

    # CDS -> gene family, keyed on coordinates (verified 363/363 join)
    cds_fam = {}
    if a.families and os.path.exists(a.families):
        for r in read_tsv(a.families):
            if int(r['family_n_loci']) >= a.min_family_loci:
                cds_fam[(r['locus_id'], r['start'], r['end'])] = r['family_id']

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

    # ---- figure geometry is settled BEFORE drawing, so the legend gets a real
    # gridspec row instead of being hung off the bottom with negative
    # coordinates (which bbox_inches='tight' then expands into a huge void) ----
    seen = {}
    for r in rows:
        if r['feature_type'] in ('IS', 'cargo'):
            seen[r['name']] = (is_colour.get(r['name']) if r['feature_type'] == 'IS'
                               else colour(r['name']), r['feature_type'])
    fam_seen = {}
    for e, f in fam_of.items():
        fam_seen.setdefault(f, set()).add(e)
    blocks = legend_blocks(fam_seen, is_colour, seen)

    PITCH = 1.75                       # vertical distance between loci
    NCOL = 6
    n_lrow = -(-len(blocks) // NCOL)
    leg_h = 0.20 * sum(max(len(b[2]) for b in blocks[i * NCOL:(i + 1) * NCOL]) + 1.45
                       for i in range(n_lrow))
    main_h = 0.95 * PITCH * len(loci) + 2.0

    fig = plt.figure(figsize=(15, main_h + leg_h))
    gs = fig.add_gridspec(2, 1, height_ratios=[main_h, leg_h],
                      hspace=0.10 + 0.9 / main_h)
    ax, lg = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    H_CDS, H_IS, H_CARGO = 0.17, 0.40, 0.30
    HEAD = a.flank / 90          # arrowhead length in bp
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
                fid = cds_fam.get((lid, r['start'], r['end']))
                ax.add_patch(Polygon(
                    gene_arrow(o, w, y, H_CDS, r['strand'], HEAD), closed=True,
                    facecolor=cds_family_colour(fid) if fid else '#f0f0f0',
                    edgecolor='#9e9e9e' if fid else '#bdbdbd',
                    lw=0.35, zorder=1))
                continue

            h = H_IS if ft == 'IS' else H_CARGO
            yy = y + 0.34 if ft == 'IS' else y - 0.36
            c = is_colour.get(nm) if ft == 'IS' else colour(nm)
            frag = r.get('cov_class') == 'fragment'
            ax.add_patch(Polygon(
                gene_arrow(o, w, yy, h, r['strand'], HEAD), closed=True,
                facecolor=c,
                edgecolor='black' if not clipped else '#555555', lw=0.9,
                linestyle='-' if not clipped else (0, (1.4, 1.0)),
                hatch='///' if frag else None, alpha=0.95, zorder=3))
            # name the vehicle and the small-but-decisive fragments
            # An element nested inside another IS is the structural event in
            # that locus, and leaving it unnamed is what made CP018961 read as
            # "ISEsa1 drawn twice in two colours": the unlabelled pale block
            # is ISKpn21, lighter only because it shares ISEsa1's family.
            if (nm.startswith('ISEsa') or (frag and ft == 'IS')
                    or (ft == 'IS' and r.get('nested_in'))):
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
                            (cx, yy + h / 2 + 0.22 + 0.30 * lvl),
                            ha='center', va='bottom', fontsize=5.6,
                            color='#222222', zorder=5, linespacing=0.95,
                            arrowprops=dict(arrowstyle='-', lw=0.4,
                                            color='#999999') if lvl else None)

        # ---- one copy, interrupted: bracket its pieces ----
        by_name = {}
        for r in feats:
            if r['feature_type'] == 'IS':
                by_name.setdefault(r['name'], []).append(r)
        interrupted = {r['nested_in'] for r in feats if r.get('nested_in')}
        for nm, parts in by_name.items():
            if nm not in interrupted or len(parts) < 2:
                continue
            parts.sort(key=lambda r: int(r['offset_from_bait']))
            x0 = (int(parts[0]['offset_from_bait'])
                  + int(parts[0]['end']) - int(parts[0]['start']) + 1)
            x1 = int(parts[-1]['offset_from_bait'])
            tot = sum(int(r['end']) - int(r['start']) + 1 for r in parts)
            yb = y + 0.34 + H_IS / 2 + 0.07
            ax.plot([x0, x0, x1, x1], [yb - 0.05, yb, yb, yb - 0.05],
                    color='#444444', lw=0.7, zorder=2)
            ax.annotate(f'one {nm}, {tot} bp in {len(parts)} pieces',
                        ((x0 + x1) / 2, yb + 0.02), ha='center', va='bottom',
                        fontsize=5.0, color='#444444', style='italic',
                        zorder=5)

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

    ytop = (len(loci) - 1) * PITCH
    for dy, lab in ((0.34, 'IS'), (0.0, 'CDS'), (-0.36, 'cargo')):
        ax.text(a.flank * 1.21, ytop + dy, lab, fontsize=6.4, va='center',
                ha='right', color='#777777', style='italic')

    ax.axvline(0, color='#1f77b4', lw=1.3, ls='--', zorder=4)
    ax.text(0, (len(loci) - 1) * PITCH + 0.95, ' relB start (bait = 0)',
            fontsize=7.5, color='#1f77b4', ha='left', va='bottom')
    # the recurring boundary, drawn so it can be checked rather than asserted
    ax.axvline(-3898, color='#ff7f0e', lw=0.9, ls=':', zorder=2)
    ax.text(-3898, -1.35, '-3898', fontsize=6.5, color='#ff7f0e', ha='center')

    ax.set_xlim(-a.flank * 1.30, a.flank * 1.22)
    ax.set_ylim(-1.55, (len(loci) - 1) * PITCH + 1.30)
    ax.set_yticks([])
    ax.set_xlabel('position relative to relB start (bp)', fontsize=9,
              labelpad=6)
    ax.set_xticks(range(-a.flank, a.flank + 1, 5000))
    ax.set_xticklabels([f'{v//1000:+d} kb' if v else '0'
                        for v in range(-a.flank, a.flank + 1, 5000)], fontsize=7.5)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.set_title(f'relBE flanking architecture, aligned on the bait'
                 f'   [palette: {a.palette}]\n'
                 'IS colour = family, shade = element within it; '
                 'cargo below each line. CDS colour = gene family; the same '
                 'colour in two rows is the same family, grey = this locus only\n'
                 'Strand is the direction each block points -- NOT whether it '
                 'sits above or below the line, which is only the track',
                 fontsize=10, pad=12)

    draw_legend(lg, blocks, NCOL)

    os.makedirs(a.outdir, exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(a.outdir, f'colour_blocks.{ext}'),
                    dpi=200, bbox_inches='tight')
    plt.close(fig)

    # ---------------- panel B: IS presence/absence heatmap ----------------
    is_rows = [r for r in rows if r['feature_type'] == 'IS']
    n_loci_with = {e: len({r['locus_id'] for r in is_rows if r['name'] == e})
                   for e in fam_of}
    # group the columns by family, families ordered by how widely they occur,
    # so the heatmap carries the same classification as the figure above
    fam_rank = {f: (-max(n_loci_with[e] for e in es), -len(es), f)
                for f, es in fam_members.items()}
    order_e = sorted(fam_of, key=lambda e: (fam_rank[fam_of[e]],
                                            -n_loci_with[e], e))

    cell = {}
    for r in is_rows:
        k = (r['locus_id'], r['name'])
        if cell.get(k) is None or r['cov_class'] == 'complete':
            cell[k] = r['cov_class']

    figb, axb = plt.subplots(figsize=(0.42 * len(order_e) + 4.2,
                                      0.46 * len(loci) + 2.6))
    CMAP = {'complete': '#08519c', 'partial': '#6baed6', 'fragment': '#c6dbef'}
    for i, lid in enumerate(loci):
        y = len(loci) - 1 - i
        for j, e in enumerate(order_e):
            axb.add_patch(Rectangle((j, y), 1, 1,
                                    facecolor=CMAP.get(cell.get((lid, e)), '#f7f7f7'),
                                    edgecolor='white', lw=1.1))
    # family band above the columns, plus a rule between families
    for f in fam_members:
        js = [j for j, e in enumerate(order_e) if fam_of[e] == f]
        j0, j1 = min(js), max(js) + 1
        axb.add_patch(Rectangle((j0 + 0.06, len(loci) + 0.12), j1 - j0 - 0.12,
                                0.34, facecolor=family_colour(f),
                                edgecolor='none', clip_on=False))
        axb.text((j0 + j1) / 2, len(loci) + 0.62, f, fontsize=6.4, ha='center',
                 va='bottom', rotation=0, color=family_colour(f),
                 fontweight='bold', clip_on=False)
        if j0:
            axb.plot([j0, j0], [0, len(loci)], color='#bbbbbb', lw=0.8)

    axb.set_xlim(0, len(order_e)); axb.set_ylim(0, len(loci))
    axb.set_xticks([j + 0.5 for j in range(len(order_e))])
    axb.set_xticklabels(order_e, rotation=90, fontsize=7)
    axb.set_yticks([len(loci) - 1 - i + 0.5 for i in range(len(loci))])
    axb.set_yticklabels([f"{l.split('__')[0]}  "
                         f"{ctx.get(l, {}).get('vehicle_class', '?')}"
                         for l in loci], fontsize=7)
    axb.set_title('IS element presence by locus, grouped by family\n'
                  f'families ordered by how widely they occur (n={len(loci)} loci)',
                  fontsize=9.5, pad=44)
    for sp in ('top', 'right', 'left', 'bottom'):
        axb.spines[sp].set_visible(False)
    axb.tick_params(length=0)
    axb.legend(handles=[Line2D([0], [0], marker='s', ls='', markersize=8,
                               markerfacecolor=v, markeredgecolor='white',
                               label=k) for k, v in CMAP.items()],
               loc='upper center', bbox_to_anchor=(0.5, -0.30), ncol=3,
               fontsize=7, frameon=False)
    figb.savefig(os.path.join(a.outdir, 'is_matrix.png'), dpi=200,
                 bbox_inches='tight')
    plt.close(figb)

    with open(os.path.join(a.outdir, 'is_matrix.tsv'), 'w') as fh:
        fh.write('locus_id\tvehicle_class\t' + '\t'.join(order_e) + '\n')
        fh.write('family\t-\t' + '\t'.join(fam_of[e] for e in order_e) + '\n')
        for lid in loci:
            fh.write(lid + '\t' + ctx.get(lid, {}).get('vehicle_class', '?'))
            for e in order_e:
                fh.write('\t' + (cell.get((lid, e)) or '0'))
            fh.write('\n')

    print(f'{len(loci)} loci, {len(is_rows)} IS calls, '
          f'{len(order_e)} elements in {len(fam_members)} families')
    for f in sorted(fam_members, key=lambda x: fam_rank[x]):
        es = sorted(fam_members[f], key=lambda e: -n_loci_with[e])
        print(f'  {f:<12} {", ".join(f"{e}({n_loci_with[e]})" for e in es)}')
    print(f'-> {a.outdir}/colour_blocks.png  .svg')
    print(f'-> {a.outdir}/is_matrix.png  is_matrix.tsv')


if __name__ == '__main__':
    main()
