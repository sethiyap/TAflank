#!/usr/bin/env python3
"""
Tier A -- bait search and sequence-level classification.

Bait-agnostic: the bait is a parameter, so this works for relBE or any other
gene of interest. Nothing here hardcodes relBE.

  bait_search.py --bait anchors/relBE_seed.fa --db blastdb/public16 \
                 --out results/bait_hits.tsv [--subregions relB:1-240,relE:240-527]

Why the filter matters rather than being a nicety: NCBI keyword searches for
relB/relE return ~188,000 records, almost all RelE/ParE superfamily background
(chromosomes carry 5-12 copies each). Sequence identity is what separates our
lineage from that background -- see PIPELINE.md F1.

PIPELINE.md rule: rejections are logged, never silent. Every HSP appears in the
output with a `tier_a` call and, where rejected, a `reason`.
"""

import argparse, os, subprocess, sys

FMT = ('qseqid sacc pident length mismatch gapopen qstart qend '
       'sstart send evalue bitscore qlen slen')
# sacc, not sseqid: -parse_seqids makes sseqid 'gb|CP018961.1|', which
# breaks every downstream join on accession.


def parse_subregions(spec):
    """'relB:1-240,relE:240-527' -> [('relB',1,240), ('relE',240,527)]"""
    out = []
    if not spec:
        return out
    for part in spec.split(','):
        name, rng = part.split(':')
        a, b = rng.split('-')
        out.append((name.strip(), int(a), int(b)))
    return out


def overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bait', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--subregions', default='')
    p.add_argument('--full-pident', type=float, default=95.0)
    p.add_argument('--full-qcov', type=float, default=95.0)
    p.add_argument('--paralogue-max-pident', type=float, default=85.0)
    p.add_argument('--evalue', default='1e-5')
    p.add_argument('--task', default='blastn',
                   help='blastn (word size 11) not megablast -- megablast '
                        'misses diverged copies at 85-92%% identity')
    p.add_argument('--blast-bin',
                   default=os.environ.get('BLAST_BIN',
                                          '/Users/pooja/opt/ncbi-blast-2.17.0+/bin'))
    a = p.parse_args()

    subs = parse_subregions(a.subregions)

    cmd = [os.path.join(a.blast_bin, 'blastn'),
           '-query', a.bait, '-db', a.db, '-task', a.task,
           '-evalue', a.evalue, '-outfmt', '6 ' + FMT,
           '-max_target_seqs', '100000', '-dust', 'no', '-soft_masking', 'false']
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        sys.exit('blastn failed:\n' + res.stderr)

    hdr = (FMT.split() + ['qcov', 'strand', 'hit_start', 'hit_end']
           + ['cov_' + s[0] for s in subs] + ['tier_a', 'reason'])
    rows, counts = [], {}

    for line in res.stdout.splitlines():
        f = line.split('\t')
        d = dict(zip(FMT.split(), f))
        pid = float(d['pident'])
        qs, qe = int(d['qstart']), int(d['qend'])
        ss, se = int(d['sstart']), int(d['send'])
        qlen = int(d['qlen'])
        qcov = round(100.0 * (qe - qs + 1) / qlen, 2)
        strand = 'minus' if ss > se else 'plus'
        lo, hi = (se, ss) if ss > se else (ss, se)

        # fraction of each named bait subregion covered by this HSP
        subcov = []
        for _, a1, a2 in subs:
            subcov.append(round(100.0 * overlap(qs, qe, a1, a2) / (a2 - a1 + 1), 1))

        reason = ''
        if pid >= a.full_pident and qcov >= a.full_qcov:
            cls = 'FULL'
        elif pid < a.paralogue_max_pident:
            cls = 'PARALOGUE'
            reason = f'pident {pid} < {a.paralogue_max_pident}: superfamily background'
        elif subs:
            # a half hit: which named subregion survived?
            present = [subs[i][0] for i, c in enumerate(subcov) if c >= 80.0]
            if len(present) == 1:
                cls = 'ORPHAN_' + present[0].upper()
                reason = f'only {present[0]} covered >=80%'
            elif not present:
                cls = 'PARTIAL'
                reason = f'qcov {qcov} below {a.full_qcov}, no subregion >=80%'
            else:
                cls = 'PARTIAL'
                reason = f'qcov {qcov} below {a.full_qcov}'
        else:
            cls = 'PARTIAL'
            reason = f'qcov {qcov} below {a.full_qcov}'

        counts[cls] = counts.get(cls, 0) + 1
        rows.append(f + [qcov, strand, lo, hi] + subcov + [cls, reason])

    with open(a.out, 'w') as fh:
        fh.write('\t'.join(hdr) + '\n')
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    print(f'{len(rows)} HSPs -> {a.out}')
    for k in sorted(counts):
        print(f'  {k:<18} {counts[k]}')


if __name__ == '__main__':
    main()
