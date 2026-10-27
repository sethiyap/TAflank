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
