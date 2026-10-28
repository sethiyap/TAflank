#!/usr/bin/env python3
"""
Import an NCBI web BLAST result and classify it exactly as a local search.

  import_blast.py --hit-table ncbi_hits.txt --bait anchors/relBE_seed.fa \
                  --out results/bait_hits.tsv \
                  --subregions relB:1-240,relE:240-527 \
                  --raw-out results/blast/01_raw.tsv \
                  --kept-out results/blast/02_tier_a.tsv

NCBI's "Hit Table (text)" and "Hit Table (csv)" downloads ARE outfmt 6 -- the
same twelve columns blastn produces by default -- so the only real work is
supplying the three the pipeline also uses:

  qlen    taken from the bait FASTA
  slen    absent; FETCH_METADATA fills it from esummary
  stitle  absent; likewise

That pairing matters: an imported run must use `--meta_source ncbi`, because
the `blast` mode parses stitle, which a hit table does not carry.

Everything downstream is unchanged. The same Tier A thresholds are applied
here as in bait_search.py, so an imported search and a local one are directly
comparable.
"""

import argparse, csv, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bait_search import (parse_subregions, overlap, write_blast6,  # noqa: E402
                         FMT_BASE)

# NCBI's field labels -> our column names
NCBI_MAP = {
    'query acc.ver': 'qseqid', 'query acc': 'qseqid', 'query id': 'qseqid',
    'subject acc.ver': 'sacc', 'subject acc': 'sacc', 'subject id': 'sacc',
    '% identity': 'pident', 'identity': 'pident',
    'alignment length': 'length', 'mismatches': 'mismatch',
    'gap opens': 'gapopen', 'q. start': 'qstart', 'q. end': 'qend',
    's. start': 'sstart', 's. end': 'send',
    'evalue': 'evalue', 'e value': 'evalue',
    'bit score': 'bitscore', 'score': 'bitscore',
    'subject length': 'slen', 'query length': 'qlen',
    'subject title': 'stitle', 'subject titles': 'stitle',
}
DEFAULT_12 = ('qseqid sacc pident length mismatch gapopen qstart qend '
              'sstart send evalue bitscore').split()


def bait_length(path):
    n = 0
    for line in open(path):
        if not line.startswith('>'):
            n += len(line.strip())
    return n


def parse_hit_table(path):
    """Return (rows, source_note). Accepts the text or csv download."""
    cols, rows, note = None, [], ''
    raw = open(path).read().splitlines()
    delim = ',' if path.lower().endswith('.csv') else '\t'
    for line in raw:
        if line.startswith('#'):
            m = re.match(r'#\s*Fields:\s*(.+)', line)
            if m:
                cols = [NCBI_MAP.get(c.strip().lower(), c.strip().lower())
                        for c in m.group(1).split(',')]
            if line.lower().startswith('# database'):
                note = line[1:].strip()
            continue
        if not line.strip():
            continue
        f = next(csv.reader([line], delimiter=delim)) if delim == ',' \
            else line.split('\t')
        rows.append(f)
    if cols is None:                       # no Fields header: assume the 12
        cols = DEFAULT_12
    return cols, rows, note


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--hit-table', required=True)
    p.add_argument('--bait', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--subregions', default='')
    p.add_argument('--raw-out', default='')
    p.add_argument('--kept-out', default='')
    p.add_argument('--full-pident', type=float, default=95.0)
    p.add_argument('--full-qcov', type=float, default=95.0)
    p.add_argument('--paralogue-max-pident', type=float, default=85.0)
    a = p.parse_args()

    subs = parse_subregions(a.subregions)
    qlen = bait_length(a.bait)
    cols, raw_rows, note = parse_hit_table(a.hit_table)
    if 'sacc' not in cols or 'pident' not in cols:
        sys.exit(f'could not recognise the columns in {a.hit_table}\n'
                 f'  parsed: {cols}\n'
                 f'  expected an NCBI "Hit Table" download (text or csv)')

    FMT = FMT_BASE
    out_cols = FMT.split()
    hdr = (out_cols + ['qcov', 'strand', 'hit_start', 'hit_end']
           + ['cov_' + s[0] for s in subs] + ['tier_a', 'reason'])
    rows, counts, raw_lines = [], {}, []

    for f in raw_rows:
        d = dict(zip(cols, f))
        d.setdefault('qlen', str(qlen))
        d.setdefault('slen', '')
        d.setdefault('stitle', '')
        # NCBI accessions carry a version; the pipeline joins on the bare form
        d['sacc'] = d.get('sacc', '').split('|')[-1]
        pid = float(d['pident'])
        qs, qe = int(d['qstart']), int(d['qend'])
        ss, se = int(d['sstart']), int(d['send'])
        qcov = round(100.0 * (qe - qs + 1) / (int(d['qlen']) or qlen), 2)
        strand = 'minus' if ss > se else 'plus'
        lo, hi = (se, ss) if ss > se else (ss, se)

        subcov = [round(100.0 * overlap(qs, qe, a1, a2) / (a2 - a1 + 1), 1)
                  for _, a1, a2 in subs]

        reason = ''
        if pid >= a.full_pident and qcov >= a.full_qcov:
            cls = 'FULL'
        elif pid < a.paralogue_max_pident:
            cls = 'PARALOGUE'
            reason = f'pident {pid} < {a.paralogue_max_pident}: superfamily background'
        elif subs:
            present = [subs[i][0] for i, c in enumerate(subcov) if c >= 80.0]
            if len(present) == 1:
                cls, reason = 'ORPHAN_' + present[0].upper(), \
                              f'only {present[0]} covered >=80%'
            else:
                cls, reason = 'PARTIAL', f'qcov {qcov} below {a.full_qcov}'
        else:
            cls, reason = 'PARTIAL', f'qcov {qcov} below {a.full_qcov}'

        counts[cls] = counts.get(cls, 0) + 1
        vals = [d.get(c, '') for c in out_cols]
        rows.append(vals + [qcov, strand, lo, hi] + subcov + [cls, reason])
        raw_lines.append('\t'.join(vals))

    with open(a.out, 'w') as fh:
        fh.write('\t'.join(hdr) + '\n')
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    if a.raw_out:
        write_blast6(a.raw_out, FMT, raw_lines,
                     'imported from NCBI web BLAST, no filtering applied',
                     f'source: {a.hit_table}  {note}')
        print(f'  [1] raw            {len(raw_lines):>7} HSPs -> {a.raw_out}')
    if a.kept_out:
        nf = len(out_cols)
        kept = ['\t'.join(str(x) for x in r[:nf])
                for r in rows if r[-2] not in ('PARALOGUE', 'REJECT')]
        write_blast6(a.kept_out, FMT, kept,
                     f'Tier A kept: pident>={a.full_pident} qcov>={a.full_qcov}; '
                     f'PARALOGUE (<{a.paralogue_max_pident}%) removed',
                     f'{len(raw_lines)} in -> {len(kept)} out')
        print(f'  [2] after Tier A   {len(kept):>7} HSPs -> {a.kept_out}')

    print(f'{len(rows)} HSPs imported -> {a.out}')
    for k in sorted(counts):
        print(f'  {k:<18} {counts[k]}')
    print('\nNOTE: a hit table carries no slen/stitle, so run the pipeline with')
    print('      --meta_source ncbi (esummary) rather than blast.')


if __name__ == '__main__':
    main()
