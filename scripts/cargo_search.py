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
    m = re.search(r'bait=(\d+)-', hdr)
    bait_start = int(m.group(1)) if m else 1

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
                     'bitscore': float(d['bitscore'])})

    # tile: best-scoring non-overlapping hits, same 30 bp tolerance as the IS layer
    hits.sort(key=lambda h: -h['bitscore'])
    kept = []
    for h in hits:
        if not any(min(h['win_end'], k['win_end'])
                   - max(h['win_start'], k['win_start']) + 1 > TOL for k in kept):
            kept.append(h)
    kept.sort(key=lambda h: h['win_start'])

    cols = ['locus_id', 'family', 'source', 'win_start', 'win_end', 'length',
            'strand', 'pident', 'cov', 'offset_from_bait', 'bitscore']
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for h in kept:
            fh.write('\t'.join(str(h[c]) for c in cols) + '\n')
    fams = sorted({h['family'] for h in kept})
    print(f'{lid}: {len(kept)} cargo hits [{",".join(fams) or "-"}] -> {a.out}')


if __name__ == '__main__':
    main()
