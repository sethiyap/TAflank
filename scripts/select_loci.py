#!/usr/bin/env python3
"""
Tier B -- join sequence hits to record metadata, compute flank availability,
and select the working set.

  select_loci.py --hits results/bait_hits.tsv --meta results/record_metadata.tsv \
                 --out results/loci_selected.tsv --flank 20000 \
                 [--tier-a FULL] [--replicon-type plasmid] [--min-flank 0]

Filters are explicit parameters, and nothing is deleted: every hit is emitted
with `selected` yes/no and a `deselect_reason`. CLAUDE.md requires that
contig-truncated loci and fragmentary hits are never silently dropped --
truncation is itself a signal, because an insertion hotspot breaks assemblies
(the blaNDM study lost 2,269 contigs exactly this way).

Flank availability mirrors scripts/relbe_window.py so the prediction here and
the sequence actually cut there cannot disagree:
  circular, replicon big enough  -> full flank both sides (window wraps)
  circular, replicon too small   -> whole replicon
  linear                         -> clamped at the ends, truncation is real
"""

import argparse


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        return hdr, [dict(zip(hdr, l.rstrip('\n').split('\t')))
                     for l in fh if l.strip()]


def flanks(hit_start, hit_end, slen, topology, flank):
    """Return (up, dn, truncated_up, truncated_dn, whole_replicon)."""
    span = hit_end - hit_start + 1
    if str(topology).lower().startswith('circ'):
        if span + 2 * flank >= slen:
            spare = slen - span
            up = spare // 2
            return up, spare - up, True, True, True
        return flank, flank, False, False, False
    up = min(flank, hit_start - 1)
    dn = min(flank, slen - hit_end)
    return up, dn, up < flank, dn < flank, False


def write_summary(path, rows, hdr, hits, a):
    """Stage-by-stage audit of the cascade.

    The filters are evaluated together, not in sequence -- a hit can fail on
    more than one -- so this reports the population at each axis and every
    distinct combination of reasons, rather than pretending to a strict
    funnel it does not implement.
    """
    ti, si, ri = hdr.index('tier_a'), hdr.index('selected'), hdr.index('replicon_type')
    def tally(i):
        c = {}
        for r in rows:
            c[r[i]] = c.get(r[i], 0) + 1
        return c
    reasons = {}
    for r in rows:
        if r[si] == 'no':
            reasons[r[-1]] = reasons.get(r[-1], 0) + 1

    out = [('input', 'blastn HSPs', len(hits), '')]
    for k, v in sorted(tally(ti).items()):
        out.append(('tier_a', k, v, 'kept' if (not a.tier_a or k in a.tier_a)
                                     else 'dropped'))
    for k, v in sorted(tally(ri).items()):
        keep = (not a.replicon_type) or k in a.replicon_type
        out.append(('replicon_type', k, v, 'kept' if keep else 'dropped'))
    for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
        out.append(('deselect_reason', k, v, 'dropped'))
    nsel = sum(1 for r in rows if r[si] == 'yes')
    out.append(('output', 'selected loci', nsel, ''))
    out.append(('params', f'tier_a={a.tier_a or "all"} '
                          f'replicon_type={a.replicon_type or "all"} '
                          f'min_flank={a.min_flank} flank={a.flank}', '', ''))

    with open(path, 'w') as fh:
        fh.write('stage\tcategory\tn\tdisposition\n')
        for st, cat, n, dis in out:
            fh.write(f'{st}\t{cat}\t{n}\t{dis}\n')
    print(f'  filter audit -> {path}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--hits', required=True)
    p.add_argument('--meta', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--flank', type=int, default=20000)
    p.add_argument('--tier-a', default='FULL',
                   help='comma-separated tier_a classes to keep')
    p.add_argument('--replicon-type', default='',
                   help='comma-separated; empty keeps all')
    p.add_argument('--min-flank', type=int, default=0,
                   help='0 keeps short-flank loci and records how much was got')
    p.add_argument('--require-complete', action='store_true')
    p.add_argument('--summary-out', default='',
                   help='write a stage-by-stage audit of the filter cascade')
    p.add_argument('--blast-out', default='',
                   help='the Tier B survivors in the same outfmt 6 columns the '
                        'earlier stages used, so the narrowing is diffable')
    a = p.parse_args()

    keep_tier = {x.strip() for x in a.tier_a.split(',') if x.strip()}
    keep_rep = {x.strip() for x in a.replicon_type.split(',') if x.strip()}

    _, hits = read_tsv(a.hits)
    _, meta = read_tsv(a.meta)
    m = {r['accession']: r for r in meta}

    hdr = ['locus_id', 'accession', 'copy', 'strand', 'tier_a',
           'pident', 'qcov', 'hit_start', 'hit_end',
           'replicon_type', 'replicon_name', 'organism', 'genus',
           'replicon_len', 'topology', 'is_complete', 'replicon_type_evidence',
           'flank_up', 'flank_dn', 'truncated_up', 'truncated_dn',
           'whole_replicon', 'window_len', 'topology_evidence',
           'selected', 'deselect_reason']
    rows = []
    percopy = {}

    for h in sorted(hits, key=lambda r: (r['sacc'], int(r['hit_start']))):
        acc = h['sacc']
        md = m.get(acc, {})
        # Unknown topology defaults to LINEAR, deliberately.
        # Treating an unknown record as circular would WRAP the window across
        # the origin, splicing the end of the sequence onto the beginning and
        # fabricating a junction that does not exist. Treating a circular one
        # as linear merely clamps and records the lost flank. One error invents
        # data; the other only loses context. BLAST output carries no topology,
        # so this is the common case in a sweep.
        topo = md.get('topology', '') or 'linear'
        topo_evid = 'record' if md.get('topology') else 'assumed_linear'
        slen = int(md.get('slen') or h.get('slen') or 0)
        hs, he = int(h['hit_start']), int(h['hit_end'])
        up, dn, tu, td, whole = flanks(hs, he, slen, topo, a.flank)

        percopy[acc] = percopy.get(acc, 0) + 1
        copy = percopy[acc]
        lid = f"{acc}__{acc}__{copy}__{h['strand']}"

        why = []
        if keep_tier and h['tier_a'] not in keep_tier:
            why.append(f"tier_a={h['tier_a']}")
        rt = md.get('replicon_type', 'unknown')
        if keep_rep and rt not in keep_rep:
            why.append(f'replicon_type={rt}')
        if a.min_flank and min(up, dn) < a.min_flank:
            why.append(f'flank {min(up, dn)} < {a.min_flank}')
        if a.require_complete and md.get('is_complete') != 'yes':
            why.append('record not complete')

        rows.append([lid, acc, copy, h['strand'], h['tier_a'],
                     h['pident'], h['qcov'], hs, he,
                     rt, md.get('replicon_name', ''), md.get('organism', ''),
                     md.get('genus', ''), slen, topo,
                     md.get('is_complete', ''),
                     md.get('replicon_type_evidence', ''),
                     up, dn, 'yes' if tu else 'no', 'yes' if td else 'no',
                     'yes' if whole else 'no',
                     (he - hs + 1) + up + dn, topo_evid,
                     'no' if why else 'yes', '; '.join(why)])

    with open(a.out, 'w') as fh:
        fh.write('\t'.join(hdr) + '\n')
        for r in rows:
            fh.write('\t'.join(str(v) for v in r) + '\n')

    sel = [r for r in rows if r[-2] == 'yes']
    nassumed = sum(1 for r in sel if r[-3] == 'assumed_linear')

    if a.summary_out:
        write_summary(a.summary_out, rows, hdr, hits, a)

    if a.blast_out:
        import os
        bcols = [c for c in ('qseqid', 'sacc', 'pident', 'length', 'mismatch',
                             'gapopen', 'qstart', 'qend', 'sstart', 'send',
                             'evalue', 'bitscore', 'qlen', 'slen', 'stitle')
                 if c in hits[0]]
        keep_acc = {(r[1], r[7], r[8]) for r in rows if r[-2] == 'yes'}
        os.makedirs(os.path.dirname(a.blast_out) or '.', exist_ok=True)
        with open(a.blast_out, 'w') as fh:
            fh.write(f'# filter: Tier B  tier_a={a.tier_a or "all"} '
                     f'replicon_type={a.replicon_type or "all"} '
                     f'min_flank={a.min_flank}\n')
            n = 0
            for h in hits:
                lo, hi = sorted((int(h['sstart']), int(h['send'])))
                if (h['sacc'], lo, hi) in keep_acc:
                    n += 1
            fh.write(f'# {len(hits)} in -> {n} out\n')
            fh.write('# ' + '\t'.join(bcols) + '\n')
            for h in hits:
                lo, hi = sorted((int(h['sstart']), int(h['send'])))
                if (h['sacc'], lo, hi) in keep_acc:
                    fh.write('\t'.join(h.get(c, '') for c in bcols) + '\n')
        print(f'  [3] after Tier B   {n:>7} HSPs -> {a.blast_out}')
    print(f'{len(rows)} hits -> {len(sel)} selected -> {a.out}')
    for r in rows:
        if r[-2] == 'no':
            print(f'  deselected {r[1]:<12} {r[-1]}')
    nt = sum(1 for r in sel if r[19] == 'yes' or r[20] == 'yes')
    if nt:
        print(f'  {nt} selected locus/loci have a short flank (kept, recorded)')
    if nassumed:
        print(f'  {nassumed} locus/loci had no topology in the metadata and were '
              f'treated as LINEAR (no wrap); flanks may understate what a '
              f'circular replicon could supply')


if __name__ == '__main__':
    main()
