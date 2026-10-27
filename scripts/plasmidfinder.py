#!/usr/bin/env python3
"""
F4 layer 5 -- replicon (Inc) typing by PlasmidFinder, via plain blastn.

  plasmidfinder.py --replicon data/public_16/fasta/CP018961.1.fna \
                   --db db/plasmidfinder/plasmidfinder.fa --out results/pf/CP018961.tsv

No CGE tooling, no conda: the PlasmidFinder database is plain FASTA and the
native BLAST already on disk types it directly. Thresholds are Acman et al.
2022 (Nat Commun 13:1131): coverage >= 80%, identity >= 95%, and overlapping
hits on one replicon resolved in favour of the higher identity.

Run on the WHOLE PARENT REPLICON, never the +/-20 kb window: Inc rep genes are
rarely within 20 kb of relBE, so a window-only run returns empty for nearly
everything.

The replicon gene is the QUERY and the replicon the SUBJECT, so coverage is
qcovs -- coverage of the rep gene. Getting this backwards is the same
qlen-vs-slen trap that once misclassified all nine loci.
"""

import argparse, os, subprocess, sys

FMT = 'qseqid qlen length pident sstart send evalue bitscore qcovs'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--replicon', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--accession', default='')
    p.add_argument('--min-cov', type=float, default=80.0)
    p.add_argument('--min-id', type=float, default=95.0)
    p.add_argument('--blast-bin',
                   default=os.environ.get('BLAST_BIN',
                                          '/Users/pooja/opt/ncbi-blast-2.17.0+/bin'))
    a = p.parse_args()

    acc = a.accession or os.path.basename(a.replicon).split('.')[0]
    cmd = [os.path.join(a.blast_bin, 'blastn'),
           '-query', a.db, '-subject', a.replicon, '-task', 'blastn',
           '-evalue', '1e-10', '-outfmt', '6 ' + FMT,
           '-dust', 'no', '-soft_masking', 'false']
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit('blastn failed:\n' + r.stderr)

    hits = []
    for line in r.stdout.splitlines():
        f = line.split('\t')
        d = dict(zip(FMT.split(), f))
        cov, pid = float(d['qcovs']), float(d['pident'])
        if cov < a.min_cov or pid < a.min_id:
            continue
        q = d['qseqid']                       # "<group>|<replicon>_<n>_<acc>"
        grp, rest = q.split('|', 1) if '|' in q else ('', q)
        ss, se = int(d['sstart']), int(d['send'])
        hits.append({'accession': acc, 'taxon_group': grp,
                     'replicon': rest.split('_')[0], 'query_id': rest,
                     'pident': pid, 'qcovs': cov, 'length': int(d['length']),
                     'start': min(ss, se), 'end': max(ss, se),
                     'strand': '-' if ss > se else '+',
                     'bitscore': float(d['bitscore'])})

    # overlapping replicon hits on the same locus -> keep the higher identity
    hits.sort(key=lambda h: (-h['pident'], -h['bitscore']))
    kept = []
    for h in hits:
        if not any(min(h['end'], k['end']) - max(h['start'], k['start']) > 0
                   for k in kept):
            kept.append(h)
    kept.sort(key=lambda h: h['start'])

    cols = ['accession', 'replicon', 'taxon_group', 'query_id', 'pident',
            'qcovs', 'length', 'start', 'end', 'strand', 'bitscore']
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for h in kept:
            fh.write('\t'.join(str(h[c]) for c in cols) + '\n')
    names = ','.join(h['replicon'] for h in kept) or '-'
    print(f'{acc}: {len(kept)} replicon(s): {names}')


if __name__ == '__main__':
    main()
