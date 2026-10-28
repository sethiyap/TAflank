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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolpath import resolve  # noqa: E402

FMT_BASE = ('qseqid sacc pident length mismatch gapopen qstart qend '
            'sstart send evalue bitscore qlen slen stitle')
# stitle comes free from any BLAST db and carries "... plasmid pX, complete
# sequence". That is enough to call plasmid vs chromosome WITHOUT esummary,
# which is what makes the search runnable on a compute node with no network
# (Gadi's standard queues). staxids/sscinames additionally need the `taxdb`
# files alongside the database -- NCBI ships them, so `update_blastdb.pl ...
# taxdb` is worth including.
FMT_TAX = ' staxids sscinames'
# sacc, not sseqid: -parse_seqids makes sseqid 'gb|CP018961.1|', which
# breaks every downstream join on accession.


def write_blast6(path, fmt, lines, applied, provenance):
    """Write BLAST outfmt 6 rows with a commented header.

    Every stage writes the SAME columns, so the files at each step of the
    cascade can be compared directly -- `wc -l` shows the narrowing and `diff`
    shows exactly which hits a filter removed.
    """
    import os
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w') as fh:
        fh.write(f'# filter: {applied}\n')
        fh.write(f'# {provenance}\n')
        fh.write('# ' + '\t'.join(fmt.split()) + '\n')
        for l in lines:
            fh.write(l.rstrip('\n') + '\n')


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
    p.add_argument('--raw-out', default='',
                   help='verbatim blastn tabular output, BEFORE any filtering')
    p.add_argument('--kept-out', default='',
                   help='the Tier A survivors, in the SAME outfmt 6 columns as '
                        '--raw-out so the two are directly comparable')
    p.add_argument('--subregions', default='')
    p.add_argument('--full-pident', type=float, default=95.0)
    p.add_argument('--full-qcov', type=float, default=95.0)
    p.add_argument('--paralogue-max-pident', type=float, default=85.0)
    p.add_argument('--evalue', default='1e-5')
    p.add_argument('--max-target-seqs', type=int, default=100000,
                   help='NOT a top-N filter -- it prunes during search, so it '
                        'must exceed the number of genuine hits expected')
    p.add_argument('--taxonomy', action='store_true',
                   help='also request staxids/sscinames; needs taxdb beside the db')
    p.add_argument('--num-threads', type=int, default=1)
    p.add_argument('--task', default='blastn',
                   help='blastn (word size 11) not megablast -- megablast '
                        'misses diverged copies at 85-92%% identity')
    p.add_argument('--blast-bin',
                   default=os.environ.get('BLAST_BIN', ''),
                   help='dir holding blastn; empty => resolve from PATH')
    a = p.parse_args()

    subs = parse_subregions(a.subregions)
    FMT = FMT_BASE + (FMT_TAX if a.taxonomy else '')

    cmd = [resolve('blastn', a.blast_bin),
           '-query', a.bait, '-db', a.db, '-task', a.task,
           '-evalue', a.evalue, '-outfmt', '6 ' + FMT,
           '-max_target_seqs', str(a.max_target_seqs),
           '-num_threads', str(a.num_threads),
           '-dust', 'no', '-soft_masking', 'false']
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        sys.exit('blastn failed:\n' + res.stderr)

    # The unprocessed search result, kept so the filtering can be audited
    # against what BLAST actually returned rather than taken on trust.
    raw_lines = res.stdout.splitlines()
    if a.raw_out:
        write_blast6(a.raw_out, FMT, raw_lines,
                     'blastn, no filtering applied',
                     ' '.join(cmd))
        print(f'  [1] raw            {len(raw_lines):>7} HSPs -> {a.raw_out}')

    hdr = ([c for c in FMT.split()] + ['qcov', 'strand', 'hit_start', 'hit_end']
           + ['cov_' + s[0] for s in subs] + ['tier_a', 'reason'])
    rows, counts = [], {}

    nfields = len(FMT.split())
    for line in res.stdout.splitlines():
        f = line.split('\t')
        if len(f) > nfields:            # stitle may itself contain a tab
            f = f[:nfields - 1] + ['\t'.join(f[nfields - 1:])]
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
        rows.append([d.get(c, '') for c in FMT.split()]
                    + [qcov, strand, lo, hi] + subcov + [cls, reason])

    with open(a.out, 'w') as fh:
        fh.write('\t'.join(hdr) + '\n')
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    # Tier A survivors, in the SAME columns as the raw dump so `diff` and
    # `wc -l` between the two show exactly what sequence identity removed.
    if a.kept_out:
        nf = len(FMT.split())
        kept = ['\t'.join(str(x) for x in r[:nf])
                for r in rows if r[-2] not in ('PARALOGUE', 'REJECT')]
        write_blast6(a.kept_out, FMT, kept,
                     f'Tier A kept: pident>={a.full_pident} qcov>={a.full_qcov} '
                     f'(FULL), orphans and partials retained; '
                     f'PARALOGUE (<{a.paralogue_max_pident}%) removed',
                     f'{len(raw_lines)} in -> {len(kept)} out')
        print(f'  [2] after Tier A   {len(kept):>7} HSPs -> {a.kept_out}')

    print(f'{len(rows)} HSPs -> {a.out}')
    for k in sorted(counts):
        print(f'  {k:<18} {counts[k]}')


if __name__ == '__main__':
    main()
