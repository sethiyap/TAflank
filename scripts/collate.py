#!/usr/bin/env python3
"""
Gather per-locus IS calls and per-record replicon types into the deliverables.

  collate.py --is-dir results/is --pf-dir results/pf --loci results/loci.tsv \
             --out-mge results/mge_hits.tsv --out-rep results/replicons.tsv \
             --out-summary results/summary.tsv

CENSUS rule: every locus in loci.tsv must appear in the summary, even if a
layer produced nothing for it, so a silent drop shows up as a visible gap
rather than a quietly smaller number.
"""

import argparse, glob, os


def read_tsv(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]


def main():
    a = argparse.ArgumentParser()
    a.add_argument('--is-dir', required=True)
    a.add_argument('--pf-dir', required=True)
    a.add_argument('--loci', required=True)
    a.add_argument('--out-mge', required=True)
    a.add_argument('--out-rep', required=True)
    a.add_argument('--out-summary', required=True)
    a = a.parse_args()

    loci = read_tsv(a.loci)

    is_rows = []
    for f in sorted(glob.glob(os.path.join(a.is_dir, '*.tsv'))):
        if f.endswith(('loci.tsv', 'summary.tsv')) or '.pf.' in f:
            continue
        is_rows += read_tsv(f)
    pf_rows = []
    for f in sorted(glob.glob(os.path.join(a.pf_dir, '*.pf.tsv'))):
        pf_rows += read_tsv(f)

    if is_rows:
        cols = list(is_rows[0].keys())
        with open(a.out_mge, 'w') as fh:
            fh.write('\t'.join(cols) + '\n')
            for r in sorted(is_rows, key=lambda x: (x['locus_id'],
                                                    int(x['win_start']))):
                fh.write('\t'.join(r.get(c, '') for c in cols) + '\n')

    if pf_rows:
        cols = list(pf_rows[0].keys())
        with open(a.out_rep, 'w') as fh:
            fh.write('\t'.join(cols) + '\n')
            for r in sorted(pf_rows, key=lambda x: (x['accession'],
                                                    int(x['start']))):
                fh.write('\t'.join(r.get(c, '') for c in cols) + '\n')

    pf_by_acc = {}
    for r in pf_rows:
        pf_by_acc.setdefault(r['accession'], []).append(r['replicon'])
    is_by_locus = {}
    for r in is_rows:
        is_by_locus.setdefault(r['locus_id'], []).append(r)

    hdr = ['locus_id', 'accession', 'organism', 'replicon_type',
           'replicon_name', 'plasmidfinder_replicons', 'replicon_len',
           'topology', 'window_len', 'flank_up', 'flank_dn',
           'truncated_up', 'truncated_dn', 'n_is', 'n_is_complete',
           'n_is_fragment', 'is_elements', 'has_ISEsa1', 'has_ISEsa2',
           'nested_events']
    out = []
    for L in loci:
        h = is_by_locus.get(L['locus_id'], [])
        els = [x['element'] for x in h]
        out.append([
            L['locus_id'], L['accession'], L['organism'], L['replicon_type'],
            L['replicon_name'],
            ','.join(sorted(set(pf_by_acc.get(L['accession'], [])))) or '-',
            L['replicon_len'], L['topology'], L['window_len'],
            L['flank_up'], L['flank_dn'], L['truncated_up'], L['truncated_dn'],
            len(h),
            sum(1 for x in h if x['cov_class'] == 'complete'),
            sum(1 for x in h if x['is_fragment'] == 'yes'),
            ','.join(els) or '-',
            'yes' if 'ISEsa1' in els else 'no',
            'yes' if 'ISEsa2' in els else 'no',
            ';'.join(f"{x['element']}>{x['nested_in']}" for x in h
                     if x.get('nested_in')) or '-',
        ])

    with open(a.out_summary, 'w') as fh:
        fh.write('\t'.join(hdr) + '\n')
        for r in out:
            fh.write('\t'.join(str(v) for v in r) + '\n')

    missing = [r[0] for r in out if r[13] == 0]
    print(f'{len(is_rows)} IS calls across {len(is_by_locus)} loci -> {a.out_mge}')
    print(f'{len(pf_rows)} replicon hits -> {a.out_rep}')
    print(f'{len(out)} loci -> {a.out_summary}')
    if missing:
        print(f'  NOTE {len(missing)} locus/loci with zero IS calls: {missing}')


if __name__ == '__main__':
    main()
