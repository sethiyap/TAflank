#!/usr/bin/env python3
"""
Merge annotation + IS calls + cargo into the figure's feature track.

  gene_blocks.py --locus loci/plasmids/X.fa --is results/is/X.is.tsv \
                 --cargo results/cargo/X.cargo.tsv --out results/blocks/X.blocks.tsv \
                 [--gff bakta/X.gff3] [--gbk data/public_16/gbk/ACC.gbk]

ANNOTATION SOURCE, in order of preference:
  1. --gff   a Bakta/Prokka GFF3 for this window          (annot_source=bakta)
  2. --gbk   the parent GenBank record, CDS clipped and
             re-coordinated into window space             (annot_source=genbank)
  3. neither -- IS and cargo only                         (annot_source=none)

Source 2 exists so this merge can be built and tested before Bakta is
available: Bakta cannot run on macOS at all (amrfinderplus has no macOS build,
any architecture), and the local loop is where iteration happens. Swapping in
Bakta later changes the source, not the logic.

TWO-TIER COLOUR KEY. A feature overlapping a tiled ISfinder call by >=50% is
coloured by ELEMENT NAME; everything else by gene family. This is not cosmetic:
annotation names IS elements by family only -- GenBank calls ISEsa1 "ISNCY
family transposase" and ISEsa2 "IS110 family transposase" -- so colouring from
annotation alone would merge the vehicle with unrelated ISNCY transposases and
make the figure misstate which element is present.
"""

import argparse, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OVERLAP_FRAC = 0.50


def read_tsv(p):
    if not p or not os.path.exists(p):
        return []
    with open(p) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]


def window_map(locus_fa):
    """Parse the locus header into what is needed to map parent coords in."""
    with open(locus_fa) as fh:
        h = fh.readline()
    g = lambda k, d='': (re.search(rf'{k}=(\S+)', h).group(1)
                         if re.search(rf'{k}=(\S+)', h) else d)
    orig = g('orig')                       # ACC:start-end  (may carry (wrap) etc)
    acc, _, span = orig.partition(':')
    nums = [int(x) for x in re.findall(r'\d+', span)]
    start, end = (nums + [0, 0])[:2]
    bait = g('bait', '1-1')
    return {'accession': acc, 'start': start, 'end': end,
            'strand': g('strand', 'plus'),
            'bait_start': int(bait.split('-')[0])}


def to_window(lo, hi, wm):
    """Parent-record coords -> window coords, honouring the reorientation."""
    if wm['strand'] == 'minus':
        return wm['end'] - hi + 1, wm['end'] - lo + 1
    return lo - wm['start'] + 1, hi - wm['start'] + 1


def cds_from_gbk(path, wm, wlen):
    from gbk_tools import parse_locus  # noqa: F401  (import check)
    txt = open(path).read()
    org = txt.find('\nORIGIN')
    feats = txt[:org] if org > 0 else txt
    out = []
    for m in re.finditer(r'^     CDS             (\S.*?)(?=^     \S|\Z)',
                         feats, re.M | re.S):
        blk = m.group(0)
        loc = m.group(1).split('\n')[0]
        nums = [int(x) for x in re.findall(r'\d+', loc)]
        if not nums:
            continue
        lo, hi = min(nums), max(nums)
        if hi < wm['start'] or lo > wm['end']:
            continue
        p = re.search(r'/product="([^"]*)"', blk, re.S)
        prod = ' '.join(p.group(1).split()) if p else 'hypothetical protein'
        g = re.search(r'/gene="([^"]*)"', blk)
        comp = 'complement' in loc
        ws, we = to_window(lo, hi, wm)
        if wm['strand'] == 'minus':
            comp = not comp
        if we < 1 or ws > wlen:
            continue
        out.append({'start': max(1, ws), 'end': min(wlen, we),
                    'strand': '-' if comp else '+',
                    'label': g.group(1) if g else '',
                    'product': prod})
    return out


def cds_from_gff(path):
    out = []
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('CDS', 'gene'):
            continue
        attr = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        out.append({'start': int(f[3]), 'end': int(f[4]), 'strand': f[6],
                    'label': attr.get('gene', ''),
                    'product': attr.get('product', 'hypothetical protein')})
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--locus', required=True)
    p.add_argument('--is', dest='is_tsv', default='')
    p.add_argument('--cargo', default='')
    p.add_argument('--gff', default='')
    p.add_argument('--gbk', default='')
    p.add_argument('--out', required=True)
    a = p.parse_args()

    lid = os.path.basename(a.locus).rsplit('.', 1)[0]
    wm = window_map(a.locus)
    wlen = sum(len(l.strip()) for l in open(a.locus) if not l.startswith('>'))

    if a.gff and os.path.exists(a.gff):
        cds, src = cds_from_gff(a.gff), 'bakta'
    elif a.gbk and os.path.exists(a.gbk):
        cds, src = cds_from_gbk(a.gbk, wm, wlen), 'genbank'
    else:
        cds, src = [], 'none'

    is_hits = read_tsv(a.is_tsv)
    cargo = read_tsv(a.cargo)

    rows = []
    # IS features first: they own their coordinates outright.
    for h in is_hits:
        s, e = int(h['win_start']), int(h['win_end'])
        rows.append([lid, 'IS', h['element'], h['family'], s, e, h['strand'],
                     h.get('cov_class', ''), h.get('is_fragment', ''),
                     h.get('clipped_at_window_edge', ''),
                     h.get('nested_in', ''), h['element'], 'isfinder', src,
                     int(h.get('offset_from_bait', s - wm['bait_start']))])

    def is_overlap(s, e):
        for h in is_hits:
            hs, he = int(h['win_start']), int(h['win_end'])
            ov = min(e, he) - max(s, hs) + 1
            if ov > 0 and ov >= OVERLAP_FRAC * (e - s + 1):
                return h['element']
        return None

    for c in cds:
        s, e = c['start'], c['end']
        el = is_overlap(s, e)
        key, srcname = (el, 'isfinder') if el else (
            (c['label'] or c['product'][:40]), 'annotation')
        rows.append([lid, 'CDS', c['label'] or c['product'][:40], '', s, e,
                     c['strand'], '', '', '', '', key, srcname, src,
                     s - wm['bait_start']])

    for c in cargo:
        s, e = int(c['win_start']), int(c['win_end'])
        rows.append([lid, 'cargo', c['family'], '', s, e, c['strand'], '', '',
                     c.get('clipped_at_window_edge', ''), '', c['family'],
                     'cargo', src, int(c['offset_from_bait'])])

    rows.sort(key=lambda r: r[4])
    cols = ['locus_id', 'feature_type', 'name', 'family', 'start', 'end',
            'strand', 'cov_class', 'is_fragment', 'clipped_at_window_edge',
            'nested_in', 'colour_key', 'colour_source', 'annot_source',
            'offset_from_bait']
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(v) for v in r) + '\n')

    # index by NAME -- r[11] was positional and pointed at colour_source only
    # until a column was inserted ahead of it
    ix = {k: i for i, k in enumerate(cols)}
    n_is = sum(1 for r in rows if r[ix['feature_type']] == 'CDS'
               and r[ix['colour_source']] == 'isfinder')
    nclip = sum(1 for r in rows
                if r[ix['clipped_at_window_edge']] not in ('', 'no'))
    note = f', {nclip} touching a window edge' if nclip else ''
    print(f'{lid}: {len(rows)} features  (annot={src}, {len(cds)} CDS, '
          f'{n_is} CDS coloured by IS element{note}) -> {a.out}')


if __name__ == '__main__':
    main()
