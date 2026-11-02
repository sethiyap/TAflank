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

# ---------------------------------------------------------------------------
# SIMPLE style. Colour stops encoding identity and encodes only ROLE, so the
# figure answers one question -- what kind of thing sits where, pointing which
# way -- instead of twenty-four. Which IS it is comes from the label, not the
# colour, so the two mustard/yellow tones carry orientation alone.
# ---------------------------------------------------------------------------
SIMPLE = {
    'is_fwd':   '#d4a017',   # mustard
    'is_rev':   '#f5dd5d',   # yellow
    'relB':     '#00b050',
    'relE':     '#e8112d',
    'czc':      '#2e75b6',   # cobalt-zinc-cadmium efflux
    'cytb':     '#7030a0',   # cytochrome b561
    'ars':      '#00b0a0',   # arsenic resistance
    'tra':      '#4a4a4a',
    'cargo':    '#c9c9c9',   # sil, intI, qnr, hok_gef: labelled, not singled out
    'cds':      '#eeeeee',
}


def simple_colour(ft, name, strand):
    if ft == 'IS':
        return SIMPLE['is_fwd'] if strand == '+' else SIMPLE['is_rev']
    return SIMPLE.get(name, SIMPLE['cargo'])

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


def consensus(rows, name, kind):
    """Median offset and length of `name` across the loci that carry it.

    Complete copies only -- a fragment's offset is where what survived starts,
    not where the element was, so averaging the two would put the consensus
    somewhere neither has ever been.
    """
    v = [(int(r['offset_from_bait']), int(r['end']) - int(r['start']) + 1,
          r['locus_id'])
         for r in rows if r['name'] == name and r['feature_type'] == kind
         and r.get('cov_class') != 'fragment']
    if not v:
        return None
    offs = sorted(o for o, _, _ in v)
    lens = sorted(l for _, l, _ in v)
    mid = len(v) // 2
    # n is LOCI, not hits -- counting distinct offsets made an element that is
    # at the same place in every locus report n=1 and lose its own emphasis
    return {'name': name, 'kind': kind, 'n': len({lid for _, _, lid in v}),
            'n_hits': len(v), 'offset': offs[mid], 'length': lens[mid],
            'lo': offs[0], 'hi': offs[-1], 'fixed': offs[0] == offs[-1]}


def draw_schematic(ax, rows, y, x0, x1, head, h, ends):
    """The idealised architecture, drawn from the data rather than asserted.

    Every block sits at the MEDIAN offset of that element across the loci
    carrying it, and an offset that is identical in every one of them is
    marked. That is the difference between a schematic and a claim: czc is
    drawn at -1184 because it is at -1184 in all four loci that have it, not
    because the slide says so.
    """
    core = [c for c in (consensus(rows, n, k) for n, k in
                        (('cytb', 'cargo'), ('czc', 'cargo'),
                         ('relB', 'cargo'), ('relE', 'cargo'),
                         ('hok_gef', 'cargo'))) if c]
    veh = [c for c in (consensus(rows, n, 'IS') for n in sorted(ends)) if c]
    tra = consensus(rows, 'tra', 'cargo')

    ax.axhline(y + h * 2.6, color='#cccccc', lw=0.8, ls='-', zorder=0)
    xw = x1 - x0
    ax.text(x0 - xw * 0.045, y, 'idealised\narchitecture', ha='right',
            va='center', fontsize=7.2, fontweight='bold', linespacing=1.25)

    # the variable region: everything upstream of the first vehicle end
    left_edge = max(x0, min([c['offset'] for c in veh] or [x0]))
    ax.add_patch(Rectangle((x0, y - h / 2), left_edge - x0, h,
                           facecolor='#f4f4f4', edgecolor='#d0d0d0',
                           lw=0.5, ls=(0, (2, 2)), zorder=1))
    ax.text((x0 + left_edge) / 2, y, 'variable — IS insertions,\n'
            'a different history in every locus', ha='center', va='center',
            fontsize=5.8, color='#777777', style='italic', linespacing=1.2)

    # cytb, czc, relB, relE and hok_gef sit inside ~2.5 kb, so their labels
    # overlap unless stacked -- downward here, away from the core bracket
    placed = []
    for c in veh + core:
        col = (SIMPLE['is_fwd'] if c['kind'] == 'IS'
               else SIMPLE.get(c['name'], SIMPLE['cargo']))
        ax.add_patch(Polygon(
            gene_arrow(c['offset'], c['length'], y, h, '+', head), closed=True,
            facecolor=col, edgecolor='black', lw=0.8, zorder=3))
        hard = c['fixed'] and c['n'] > 2
        txt = (f"{c['name']}  {c['offset']:+d}"
               + ('  (same in all ' + str(c['n']) + ')' if hard else ''))
        # reserve width from the TEXT, not a fixed span -- "(same in all 4)"
        # is three times as wide as a bare name and overran its neighbour
        cx, half = c['offset'] + c['length'] / 2, len(txt) * 150 / 2
        lvl = 0
        while any(abs(cx - px) < half + ph and lv == lvl
                  for px, ph, lv in placed):
            lvl += 1
        placed.append((cx, half, lvl))
        ax.annotate(txt,
                    (cx, y - h / 2 - 0.08 - 0.26 * lvl),
                    ha='center', va='top', fontsize=5.4,
                    color='#111111' if hard else '#555555',
                    fontweight='bold' if hard else 'normal',
                    arrowprops=dict(arrowstyle='-', lw=0.4,
                                    color='#bbbbbb') if lvl else None)
    depth = max([lv for _, _, lv in placed] or [0])

    ax.text(x0 - xw * 0.045, y - 0.30 - 0.26 * depth,
            'median offset across the loci carrying each element',
            ha='left', va='top', fontsize=5.2, color='#999999', style='italic')

    if tra and tra['lo'] < x1:
        ax.add_patch(Rectangle((tra['lo'], y - h / 2),
                               min(x1, tra['hi'] + tra['length']) - tra['lo'],
                               h, facecolor=SIMPLE['tra'], alpha=0.30,
                               edgecolor='#4a4a4a', lw=0.5, ls=(0, (2, 2)),
                               zorder=1))
        ax.text((tra['lo'] + min(x1, tra['hi'] + tra['length'])) / 2, y,
                f"tra — conjugation machinery\n{tra['lo']:+d} to {tra['hi']:+d}",
                ha='center', va='center', fontsize=5.8, color='#444444',
                style='italic', linespacing=1.2)

    # the conserved core bracket
    if core:
        x0 = min(c['offset'] for c in core)
        x1 = max(c['offset'] + c['length'] for c in core)
        yb = y + h / 2 + 0.14
        ax.plot([x0, x0, x1, x1], [yb - 0.06, yb, yb, yb - 0.06],
                color='#111111', lw=0.9, zorder=4)
        ax.annotate('conserved core — the extended anchor '
                    '(CLAUDE.md Stage 2)', ((x0 + x1) / 2, yb + 0.03),
                    ha='center', va='bottom', fontsize=6.0,
                    fontweight='bold', color='#111111')
    return depth


def label_plan(feats, which, ends, min_span=1300, kinds=('IS',),
               skip=()):
    """Which IS blocks get a label, and on which stacked level.

    One rule for every hit instead of three special cases. The old condition
    was `name.startswith('ISEsa') or fragment or nested`, which hardcoded this
    project's vehicle: run the pipeline with a different bait and the elements
    that mattered would silently go unlabelled, exactly the failure that
    `--vehicle_ends` was added to classify_context.py to prevent.

    `which`:
      all      every IS -- an unlabelled block is what made ISKpn21 read as
               a second colour of ISEsa1
      distinct one label per distinct element (complete and fragment counted
               apart). The default in the simple style, where colour no
               longer tells elements apart so repeats add nothing
      key     the vehicle ends, fragments, and anything nested in another IS
      none    no labels

    `skip` drops names whose colour already identifies them uniquely. relB and
    relE are green and red and sit at the same place in every row; labelling
    them bold forced two extra stack levels at x=0 on all 14 rows, which set
    the row pitch for the whole figure.

    Returns {index in feats: level}. Level is the greedy stack height needed
    to clear every label already placed on that row, so the caller can size
    the row before drawing into it.
    """
    out, placed, seen_name = {}, [], set()
    if which == 'none':
        return out
    for i, r in enumerate(feats):
        if r['feature_type'] not in kinds or r['name'] in skip:
            continue
        frag = r.get('cov_class') == 'fragment'
        if which == 'key' and not (r['name'] in ends or frag
                                   or r.get('nested_in')):
            continue
        if which == 'distinct':
            # one label per distinct element, but a fragment is kept separate
            # from a complete copy of the same element -- in CP018961 that is
            # the difference between the 1540 bp ISEsa1 and the 64 bp remnant,
            # which is the evidence, not a repeat
            key = (r['name'], frag)
            if key in seen_name:
                continue
            seen_name.add(key)
        o = int(r['offset_from_bait'])
        w = int(r['end']) - int(r['start']) + 1
        cx, half = o + w / 2, max(w, min_span) / 2
        lvl = 0
        while any(abs(cx - px) < (half + ph) and lv == lvl
                  for px, ph, lv in placed):
            lvl += 1
        placed.append((cx, half, lvl))
        out[i] = lvl
    return out


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
    p.add_argument('--no-label-names', default='relB,relE',
                   help='features never labelled on the figure because their '
                        'colour is already unique. Empty string labels all.')
    p.add_argument('--xlim', default='',
                   help='display range as LOW,HIGH in bp, e.g. -8000,3000. '
                        'Zooms the FIGURE only -- the windows stay whatever '
                        'they were extracted at, so nothing is re-analysed '
                        'and no call changes. Features outside are counted '
                        'in the title, never silently hidden.')
    p.add_argument('--no-schematic', action='store_true',
                   help='omit the idealised-architecture panel at the bottom')
    p.add_argument('--style', choices=['simple', 'tracks'], default='simple',
                   help="'simple' = one row per locus, colour by role, IS "
                        "mustard forward / yellow reverse; 'tracks' = IS, CDS "
                        "and cargo on separate lines, colour by family")
    p.add_argument('--label-is',
                   choices=['all', 'distinct', 'key', 'none'], default='',
                   help="which IS blocks get a name drawn on the figure")
    p.add_argument('--vehicle-ends', default='ISEsa1,ISEsa2',
                   help='the elements bracketing the core -- emphasised in '
                        'the figure. Must match ENDS / --vehicle_ends used '
                        'for the run, or the figure highlights the wrong ones')
    p.add_argument('--palette', choices=sorted(PALETTES), default='deck',
                   help="'deck' = relBE_schematic.pptx slide 11 (default), "
                        "'lab2018' = the Sept 2018 lab-meeting deck, "
                        "'default' = this script's own palette")
    p.add_argument('--min-family-loci', type=int, default=2,
                   help='colour a gene family only if it spans this many loci '
                        '(1 would give every singleton its own colour and turn '
                        'the track into confetti)')
    p.add_argument('--groups', default='',
                   help='results/structural_groups.tsv (structural_groups.py) '
                        '-- paints a swatch beside each locus label; loci '
                        'with matching swatches share >=90%% core '
                        'architecture and plausibly share one history')
    a = p.parse_args()

    simple = a.style == 'simple'
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
    groups = {}
    if a.groups and os.path.exists(a.groups):
        for r in read_tsv(a.groups):
            groups[r['locus_id']] = r

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
    if simple:
        blocks = [
            ('insertion sequences', '#333333', [
                ('IS, forward', dict(facecolor=SIMPLE['is_fwd'])),
                ('IS, reverse', dict(facecolor=SIMPLE['is_rev'])),
                ('fragment <50% cov', dict(facecolor=SIMPLE['is_fwd'],
                                           hatch='///')),
                ('runs off the window', dict(facecolor=SIMPLE['is_fwd'],
                                             linestyle=(0, (1.4, 1.0)))),
            ]),
            ('the TA pair', '#333333', [
                ('relB', dict(facecolor=SIMPLE['relB'])),
                ('relE', dict(facecolor=SIMPLE['relE'])),
            ]),
            ('the conserved core', '#333333', [
                ('cytb (cytochrome b561)', dict(facecolor=SIMPLE['cytb'])),
                ('czc (Co-Zn-Cd efflux)', dict(facecolor=SIMPLE['czc'])),
            ]),
            ('other context', '#333333', [
                ('ars (arsenic res.)', dict(facecolor=SIMPLE['ars'])),
                ('tra (conjugation)', dict(facecolor=SIMPLE['tra'])),
                ('sil, intI, qnr, hok_gef', dict(facecolor=SIMPLE['cargo'])),
                ('CDS', dict(facecolor=SIMPLE['cds'], edgecolor='#cccccc')),
            ]),
        ]
    else:
        blocks = legend_blocks(fam_seen, is_colour, seen)

    # Structural groups get their OWN full-width row rather than a slot in the
    # multi-column swatch grid: a group's member list ("6 loci: CP012599,
    # CP018961, ...") is long text that the narrow per-column layout above
    # was built for short element names, not for.
    group_rows = []
    if groups:
        seen_gid = set()
        for lid in loci:
            g = groups.get(lid)
            if not g or g['group_colour'] == '#cccccc' or g['group_id'] in seen_gid:
                continue
            seen_gid.add(g['group_id'])
            members = [m.split('__')[0] for m in g['group_members'].split(',')]
            group_rows.append((g['group_colour'], len(members),
                               ', '.join(members), g['shared_core']))

    if not a.label_is:
        a.label_is = 'distinct' if simple else 'all'
    ends = {e.strip() for e in a.vehicle_ends.split(',') if e.strip()}
    kinds = ('IS', 'cargo') if simple else ('IS',)
    skip = {n.strip() for n in a.no_label_names.split(',') if n.strip()}
    plans = {l: label_plan(by_locus[l], a.label_is, ends, kinds=kinds,
                           skip=skip)
             for l in loci}
    deepest = max([lv for p in plans.values() for lv in p.values()] or [0])
    # rows must be tall enough for the deepest label stack, or labels from one
    # locus climb into the row above
    # Rows sit as close as the deepest label stack allows. The +1 is the
    # split-copy bracket, which is drawn above the labels; without a bracket
    # anywhere that headroom is wasted, so only pay for it when one exists.
    has_bracket = any(
        r.get('nested_in') for rs in by_locus.values() for r in rs)
    PITCH = ((0.52 if simple else 1.20)
             + 0.26 * (deepest + (1 if has_bracket else 0)))
    NCOL = 6
    n_lrow = -(-len(blocks) // NCOL)
    leg_h = 0.20 * sum(max(len(b[2]) for b in blocks[i * NCOL:(i + 1) * NCOL]) + 1.45
                       for i in range(n_lrow))
    main_h = 0.78 * PITCH * (len(loci) + (0 if a.no_schematic else 1.7)) + 2.0
    grp_h = 0.34 * (len(group_rows) + 1.3) if group_rows else 0

    fig = plt.figure(figsize=(15, main_h + leg_h + grp_h))
    if group_rows:
        gs = fig.add_gridspec(3, 1, height_ratios=[main_h, grp_h, leg_h],
                              hspace=0.08 + 0.7 / main_h)
        ax, gp, lg = (fig.add_subplot(gs[0]), fig.add_subplot(gs[1]),
                     fig.add_subplot(gs[2]))
    else:
        gs = fig.add_gridspec(2, 1, height_ratios=[main_h, leg_h],
                          hspace=0.10 + 0.9 / main_h)
        ax, lg = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        gp = None

    # simple: one row -- CDS a thin underlay, IS and cargo full height on it
    H_CDS, H_IS, H_CARGO = (0.12, 0.40, 0.40) if simple else (0.17, 0.40, 0.30)
    IS_DY, CARGO_DY = (0.0, 0.0) if simple else (0.34, -0.36)
    HEAD = a.flank / 90          # arrowhead length in bp
    # Every label anchors to the DISPLAY range, not the extraction flank.
    # Anchoring to the flank put the locus names 20 kb off-screen in a zoomed
    # figure, and bbox_inches='tight' then grew the canvas to include them --
    # a 6161 px image of an 11 kb window.
    if a.xlim:
        X0, X1 = (int(v) for v in a.xlim.split(','))
    else:
        X0, X1 = -a.flank, a.flank
    XW = X1 - X0
    label_rows = []   # (lid, y, text object) -- swatches are placed after
                       # axis limits are final, see the note below the loop
    for i, lid in enumerate(loci):
        y = (len(loci) - 1 - i) * PITCH
        feats = by_locus[lid]
        plan = plans[lid]
        offs = [int(r['offset_from_bait']) for r in feats]
        lo, hi = max(X0, min(offs + [X0])), min(X1, max(offs + [X1]))
        ax.plot([lo, hi], [y, y], color='#d9d9d9', lw=1.0, zorder=0,
                solid_capstyle='butt')

        for fi, r in enumerate(feats):
            o = int(r['offset_from_bait'])
            w = int(r['end']) - int(r['start']) + 1
            # Clip to the display range. Without this a feature at -12 kb is
            # still drawn in a figure zoomed to -8 kb, and matplotlib parks it
            # against the axes edge -- on top of the locus labels, looking
            # like data that is actually off-screen.
            vo, ve = max(o, X0), min(o + w, X1)
            if ve <= vo:
                plan.pop(fi, None)
                continue
            o, w = vo, ve - vo
            ft, nm = r['feature_type'], r['name']
            clipped = r.get('clipped_at_window_edge', 'no') not in ('no', '-', '')

            if ft == 'CDS':
                fid = cds_fam.get((lid, r['start'], r['end']))
                ax.add_patch(Polygon(
                    gene_arrow(o, w, y, H_CDS, r['strand'], HEAD), closed=True,
                    facecolor=(SIMPLE['cds'] if simple else
                               cds_family_colour(fid) if fid else '#f0f0f0'),
                    edgecolor='#cccccc' if simple else
                              ('#9e9e9e' if fid else '#bdbdbd'),
                    lw=0.3 if simple else 0.35, zorder=1))
                continue

            h = H_IS if ft == 'IS' else H_CARGO
            yy = y + (IS_DY if ft == 'IS' else CARGO_DY)
            c = (simple_colour(ft, nm, r['strand']) if simple else
                 is_colour.get(nm) if ft == 'IS' else colour(nm))
            frag = r.get('cov_class') == 'fragment'
            ax.add_patch(Polygon(
                gene_arrow(o, w, yy, h, r['strand'], HEAD), closed=True,
                facecolor=c,
                edgecolor='black' if not clipped else '#555555', lw=0.9,
                linestyle='-' if not clipped else (0, (1.4, 1.0)),
                hatch='///' if frag else None, alpha=0.95, zorder=3))
            # name the vehicle and the small-but-decisive fragments
            if fi in plan:
                lvl = plan[fi]
                emph = nm in ends or nm in ('relB', 'relE')
                ax.annotate(
                    f'{nm}\n{w} bp' if frag else nm,
                    (o + w / 2, y + (IS_DY if ft == 'IS' else CARGO_DY)
                     + h / 2 + 0.12 + 0.26 * lvl),
                    ha='center', va='bottom',
                    fontsize=6.0 if emph else 5.4,
                    fontweight='bold' if emph else 'normal',
                    color='#111111' if emph else '#444444',
                    zorder=5, linespacing=0.95,
                    arrowprops=dict(arrowstyle='-', lw=0.4,
                                    color='#aaaaaa') if lvl else None)

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
            # above the deepest label on this row, or it lands on top of them
            yb = (y + IS_DY + H_IS / 2 + 0.12
                  + 0.26 * (max(plan.values(), default=0) + 1))
            ax.plot([x0, x0, x1, x1], [yb - 0.07, yb, yb, yb - 0.07],
                    color='#444444', lw=0.7, zorder=2)
            ax.annotate(f'one {nm}, {tot} bp in {len(parts)} pieces',
                        ((x0 + x1) / 2, yb + 0.02), ha='center', va='bottom',
                        fontsize=5.0, color='#444444', style='italic',
                        zorder=5)

        if not any(r['feature_type'] == 'CDS' for r in feats):
            ax.text(X0 - XW * 0.02, y - 0.62, 'no CDS in the GenBank record',
                    fontsize=5.8, style='italic', color='#aa0000',
                    ha='left', va='center')
        cls = ctx.get(lid, {}).get('vehicle_class', '?')
        acc = lid.split('__')[0]
        org = meta.get(acc, {}).get('organism', '')
        org = ' '.join(org.split()[:2]) if org else ''
        label_obj = ax.text(X0 - XW * 0.045, y, f'{acc}\n{org}', ha='right',
                            va='center', fontsize=7.2, linespacing=1.25,
                            style='italic' if org else 'normal')
        label_rows.append((lid, y, label_obj))
        ax.text(X1 + XW * 0.02, y, cls, ha='left', va='center', fontsize=7,
                color={'COMPLETE': '#1a7f1a', 'DECAYED': '#b8860b',
                       'HALF-CORE': '#d2691e'}.get(cls, '#999999'),
                fontweight='bold')

    ytop = (len(loci) - 1) * PITCH
    if not simple:
        for dy, lab in ((0.34, 'IS'), (0.0, 'CDS'), (-0.36, 'cargo')):
            ax.text(X1 + XW * 0.21, ytop + dy, lab, fontsize=6.4, va='center',
                    ha='right', color='#777777', style='italic')

    y_sch = -1.55 * PITCH
    sch_depth = 0
    if not a.no_schematic:
        sch_depth = draw_schematic(ax, rows, y_sch, X0, X1, HEAD, 0.40, ends)

    ax.axvline(0, color='#1f77b4', lw=1.3, ls='--', zorder=4)
    ax.text(0, (len(loci) - 1) * PITCH + 0.30 + 0.26 * (deepest + 1),
            ' relB start (bait = 0)', fontsize=7.5, color='#1f77b4',
            ha='left', va='bottom')
    # the recurring boundary, drawn so it can be checked rather than asserted
    ax.axvline(-3898, color='#ff7f0e', lw=0.9, ls=':', zorder=2)
    ybot = (y_sch - 0.95 - 0.26 * sch_depth) if not a.no_schematic else -1.35
    ax.text(-3898, ybot, '-3898', fontsize=6.5, color='#ff7f0e', ha='center')

    zoom_note = ''
    if a.xlim:
        lo_x, hi_x = (int(v) for v in a.xlim.split(','))
        span = hi_x - lo_x
        # Count what the zoom puts off-screen and say so. A figure that
        # silently drops a third of its IS calls looks like a finding.
        n_out = sum(1 for r in rows
                    if r['feature_type'] in ('IS', 'cargo')
                    and (int(r['offset_from_bait'])
                         + int(r['end']) - int(r['start']) + 1 < lo_x
                         or int(r['offset_from_bait']) > hi_x))
        ax.set_xlim(lo_x - span * 0.32, hi_x + span * 0.14)
        step = 1000 if span <= 12000 else 5000
        ticks = list(range(-(-lo_x // step) * step, hi_x + 1, step))
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{v/1000:+.0f} kb' if v else '0' for v in ticks],
                           fontsize=7.5)
        zoom_note = (f'   zoom {lo_x:+d}..{hi_x:+d} of ±{a.flank // 1000} kb'
                     f' — {n_out} IS/cargo outside')
    else:
        ax.set_xlim(-a.flank * 1.30, a.flank * 1.22)
    ax.set_ylim(ybot - 0.45,
                (len(loci) - 1) * PITCH + 0.62 + 0.26 * (deepest + 1))

    # Group swatches are placed AFTER xlim/ylim are final and by MEASURING the
    # label text's own rendered width, not by guessing a gap as a fraction of
    # the x-range. A fraction of XW is a fixed physical gap only at the flank
    # this was tuned against; --xlim shrinks XW but the font size (and so the
    # text's rendered width) does not shrink with it, so that gap collapsed
    # to less than the text width and the swatch was drawn ON TOP of
    # "Escherichia coli" rather than beside it.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for lid, y, label_obj in label_rows:
        grp = groups.get(lid)
        if not grp or grp['group_colour'] == '#cccccc':
            continue
        bbox = label_obj.get_window_extent(renderer=renderer)
        y_px = ax.transData.transform((0, y))[1]
        inv = ax.transData.inverted()
        (x_right, _) = inv.transform((bbox.x0 - 5, y_px))   # 5 px gap
        (x_left, _) = inv.transform((bbox.x0 - 15, y_px))   # 10 px wide
        ax.add_patch(Rectangle((x_left, y - 0.17), x_right - x_left, 0.34,
                               facecolor=grp['group_colour'], edgecolor='none',
                               clip_on=False))
    ax.set_yticks([])
    ax.set_xlabel('position relative to relB start (bp)', fontsize=9,
              labelpad=6)
    if not a.xlim:
        ax.set_xticks(range(-a.flank, a.flank + 1, 5000))
        ax.set_xticklabels([f'{v//1000:+d} kb' if v else '0'
                            for v in range(-a.flank, a.flank + 1, 5000)],
                           fontsize=7.5)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    # A figure showing 7 of 14 loci without saying so overstates its scope.
    # If every locus drawn shares one replicon type, that is a filter, not a
    # finding, and the title has to say which.
    kinds = {ctx.get(l, {}).get('replicon_type', '') for l in loci}
    # (scope is built here and may be extended by --xlim below)
    scope = (f'   [{kinds.pop()}s only — {len(loci)} loci]'
             if len(kinds) == 1 and next(iter(kinds), '') else
             f'   [{len(loci)} loci]')
    ax.set_title(('relBE flanking architecture, aligned on the bait'
                  + scope + zoom_note + '\n'
                  'IS mustard = forward, yellow = reverse (type is in the '
                  'label, not the colour); strand is the way each block points'
                  if simple else
                  f'relBE flanking architecture, aligned on the bait'
                  f'   [palette: {a.palette}]\n'
                 'IS colour = family, shade = element within it; '
                 'cargo below each line. CDS colour = gene family; the same '
                 'colour in two rows is the same family, grey = this locus only\n'
                 'Strand is the direction each block points -- NOT whether it '
                 'sits above or below the line, which is only the track'),
                 fontsize=10, pad=12)

    if gp is not None:
        gp.set_xlim(0, 1); gp.set_ylim(len(group_rows) + 1, 0); gp.axis('off')
        gp.text(0, 0.3, '≥ 90% core architecture (same vehicle + '
                'context genes) — plausibly one history per colour',
                fontsize=7.6, fontweight='bold', va='center')
        for i, (col, n, members, shared) in enumerate(group_rows):
            yy = i + 1.15
            gp.add_patch(Rectangle((0.0, yy - 0.28), 0.014, 0.52,
                                   facecolor=col, edgecolor='none'))
            gp.text(0.022, yy, f'{n} loci:  {members}', fontsize=7,
                    va='center')
            gp.text(0.70, yy, f'shared: {shared}', fontsize=6.2,
                    va='center', color='#666666')

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
