#!/usr/bin/env python3
"""
Tier C input -- which cargo genes sit in a locus window, and at what offset.

  cargo_search.py --locus loci/plasmids/X.fa --db assets/cargo_genes.fa \
                  --out results/cargo/X.cargo.tsv

Deliberately NOT an annotation pass. Characterising context needs only the
handful of gene families the biology turns on, and blastn against a 333-sequence
curated set answers that in milliseconds. Bakta is for the figure's complete
gene track and runs later, on far fewer loci.

Cargo presence is RECORDED, never used to admit or exclude a locus -- see
PIPELINE_DESCRIPTION.md section 4b.

Same window-edge caveat as the IS layer, for the same reason: a gene lying
across the boundary of the +/-20 kb window reports low `cov` because we stopped
looking, not because the gene is broken. The edge columns separate the two
causes -- see is_search.py's docstring for the full argument.
"""

import argparse, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolpath import resolve  # noqa: E402

FMT = 'qseqid pident length qlen qstart qend sstart send evalue bitscore'
TOL = 30


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--locus', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--min-pident', type=float, default=80.0)
    p.add_argument('--min-cov', type=float, default=50.0)
    p.add_argument('--blast-bin', default=os.environ.get('BLAST_BIN', ''))
    a = p.parse_args()

    lid = os.path.basename(a.locus).rsplit('.', 1)[0]
    with open(a.locus) as fh:
        hdr = fh.readline()
        win_len = sum(len(l.strip()) for l in fh)
    m = re.search(r'bait=(\d+)-', hdr)
    bait_start = int(m.group(1)) if m else 1
    # post-reorientation, so truncated_up is the LEFT edge of this FASTA
    trunc_left = 'truncated_up=yes' in hdr
    trunc_right = 'truncated_dn=yes' in hdr

    r = subprocess.run([resolve('blastn', a.blast_bin), '-query', a.db,
                        '-subject', a.locus, '-task', 'blastn',
                        '-evalue', '1e-10', '-outfmt', '6 ' + FMT,
                        '-dust', 'no', '-soft_masking', 'false'],
                       capture_output=True, text=True)
    if r.returncode:
        sys.exit('blastn failed:\n' + r.stderr)

    hits = []
    for line in r.stdout.splitlines():
        d = dict(zip(FMT.split(), line.split('\t')))
        pid, L, qlen = float(d['pident']), int(d['length']), int(d['qlen'])
        cov = 100.0 * L / qlen
        if pid < a.min_pident or cov < a.min_cov:
            continue
        fam, acc, lt, coords = (d['qseqid'].split('|') + ['', '', ''])[:4]
        ss, se = int(d['sstart']), int(d['send'])
        lo, hi = min(ss, se), max(ss, se)
        hits.append({'locus_id': lid, 'family': fam, 'source': d['qseqid'],
                     'win_start': lo, 'win_end': hi,
                     'strand': '-' if ss > se else '+',
                     'pident': round(pid, 2), 'cov': round(cov, 1),
                     'length': L, 'offset_from_bait': lo - bait_start,
                     'qlen': qlen, 'qstart': int(d['qstart']),
                     'qend': int(d['qend']),
                     'bitscore': float(d['bitscore'])})

    # tile: best-scoring non-overlapping hits, same 30 bp tolerance as the IS layer
    hits.sort(key=lambda h: -h['bitscore'])
    kept = []
    for h in hits:
        if not any(min(h['win_end'], k['win_end'])
                   - max(h['win_start'], k['win_start']) + 1 > TOL for k in kept):
            kept.append(h)
    kept.sort(key=lambda h: h['win_start'])

    # ---- window-edge clipping (see module docstring) ----
    for h in kept:
        qs, qe = sorted((h['qstart'], h['qend']))
        if h['strand'] == '+':
            miss_l, miss_r = qs - 1, h['qlen'] - qe
        else:
            miss_l, miss_r = h['qlen'] - qe, qs - 1
        sides, causes, missing = [], [], 0
        if h['win_start'] <= 1 and miss_l > 0:
            sides.append('left')
            causes.append('replicon' if trunc_left else 'window')
            missing += miss_l
        if h['win_end'] >= win_len and miss_r > 0:
            sides.append('right')
            causes.append('replicon' if trunc_right else 'window')
            missing += miss_r
        h['clipped_at_window_edge'] = ('both' if len(sides) == 2
                                       else sides[0] if sides else 'no')
        h['clip_cause'] = ('-' if not causes else
                           causes[0] if len(set(causes)) == 1 else 'mixed')
        h['clipped_missing_bp'] = missing

    cols = ['locus_id', 'family', 'source', 'win_start', 'win_end', 'length',
            'strand', 'pident', 'cov', 'clipped_at_window_edge', 'clip_cause',
            'clipped_missing_bp', 'offset_from_bait', 'bitscore']
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for h in kept:
            fh.write('\t'.join(str(h[c]) for c in cols) + '\n')
    fams = sorted({h['family'] for h in kept})
    nclip = sum(1 for h in kept if h['clip_cause'] == 'window')
    note = f', {nclip} cut by the window edge' if nclip else ''
    print(f'{lid}: {len(kept)} cargo hits '
          f'[{",".join(fams) or "-"}]{note} -> {a.out}')


if __name__ == '__main__':
    main()
