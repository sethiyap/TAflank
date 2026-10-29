#!/usr/bin/env python3
"""
Tier C -- classify each locus by what is around it.

  classify_context.py --loci results/loci.tsv --is-dir results/is \
                      --cargo-dir results/cargo --out results/context.tsv

TWO INDEPENDENT AXES, deliberately never merged (PIPELINE_DESCRIPTION.md 4b):

  lineage  -- from SEQUENCE identity alone, decided upstream in Tier A.
              A paralogue is a DIFFERENT GENE, not a copy without a vehicle.

  context  -- from what is found in the window. COMPLETE / DECAYED /
              HALF-CORE / NO CORE for the vehicle, plus a cargo profile.

Collapsing them -- "no vehicle, therefore a paralogue" -- would misclassify
CP017220 and CP006642, both 99.4% identical to the seed and both genuinely the
lineage, and would delete the evidence for repeated independent acquisition.

Nothing here filters. Every locus gets a row, including loci with no IS and no
cargo at all.
"""

import argparse, glob, os, sys


def read_tsv(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]


def vehicle_class(is_rows, ends=('ISEsa1', 'ISEsa2')):
    """COMPLETE / DECAYED / HALF-CORE / NO CORE, following PIPELINE.md F3.

    Returns (class, ends_present, interrupted, edge_clipped).
    """
    present, interrupted, clipped = {}, set(), {}
    for h in is_rows:
        e = h['element']
        if e not in ends:
            continue
        cov_h = float(h.get('query_cov') or 0)
        if cov_h >= clipped.get(e, (-1.0, False))[0]:
            clipped[e] = (cov_h, h.get('clip_cause') == 'window')
        # BEST coverage per element, not per HSP. An element split in two by a
        # later insertion contributes a large piece and a small one; judging on
        # the small one alone would call a 96%-covered ISEsa1 "decayed" purely
        # because its 64 bp tail is also reported.
        present[e] = max(present.get(e, 0.0), float(h.get('query_cov') or 0))
    for h in is_rows:
        if h.get('nested_in') in ends:
            interrupted.add(h['nested_in'])

    # An end whose best HSP was cut by the window edge has an UNKNOWN true
    # coverage -- the rest of the element lies outside what we extracted. Low
    # coverage there is not decay, so it must not drive the DECAYED call. The
    # element is reported as edge_clipped instead: a prompt to widen the window,
    # not a conclusion about the vehicle.
    edge_clipped = {e for e, (_, c) in clipped.items() if c}
    degraded = {e for e, cov in present.items()
                if cov < 90 and e not in edge_clipped}

    ec = ','.join(sorted(edge_clipped))
    n = len(present)
    if n == 0:
        return 'NO CORE', '', '', ec
    if n == 1:
        return ('HALF-CORE', next(iter(present)),
                ';'.join(sorted(interrupted)), ec)
    cls = 'DECAYED' if degraded else 'COMPLETE'
    return cls, ','.join(sorted(present)), ';'.join(sorted(interrupted)), ec


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--loci', required=True)
    p.add_argument('--is-dir', required=True)
    p.add_argument('--cargo-dir', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ends', default='ISEsa1,ISEsa2',
                   help='comma-separated IS elements that bracket the core. '
                        'These are lineage-specific: ISEsa1/ISEsa2 belong to '
                        'the relBE metal-resistance lineage, so a DIFFERENT '
                        'bait almost certainly needs different ends, or every '
                        'locus comes back NO CORE.')
    a = p.parse_args()

    loci = read_tsv(a.loci)
    is_by, cargo_by = {}, {}
    ends = tuple(x.strip() for x in a.ends.split(',') if x.strip())
    if not ends:
        sys.exit('--ends must name at least one element')

    for f in glob.glob(os.path.join(a.is_dir, '*.is.tsv')):
        for r in read_tsv(f):
            is_by.setdefault(r['locus_id'], []).append(r)
    for f in glob.glob(os.path.join(a.cargo_dir, '*.cargo.tsv')):
        for r in read_tsv(f):
            cargo_by.setdefault(r['locus_id'], []).append(r)

    hdr = ['locus_id', 'accession', 'replicon_type', 'vehicle_class',
           'vehicle_ends', 'vehicle_interrupted', 'vehicle_edge_clipped',
           'n_is', 'n_is_fragment',
           'n_is_fragment_real', 'cargo_families',
           'has_czc', 'has_cytb', 'has_hok_gef',
           'czc_offset', 'hok_gef_offset', 'metal_module', 'note']
    rows = []
    for L in loci:
        lid = L['locus_id']
        ish = is_by.get(lid, [])
        car = cargo_by.get(lid, [])
        vcls, vends, vint, vclip = vehicle_class(ish, ends=ends)
        fams = sorted({c['family'] for c in car})
        off = {}
        for c in car:
            off.setdefault(c['family'], c['offset_from_bait'])
        czc, cytb = 'czc' in fams, 'cytb' in fams
        rows.append([
            lid, L['accession'], L.get('replicon_type', ''),
            vcls, vends, vint, vclip or '-',
            len(ish), sum(1 for h in ish if h.get('is_fragment') == 'yes'),
            # fragments the window edge cut are excluded: they say nothing
            # about decay, only about where we stopped looking
            sum(1 for h in ish if h.get('is_fragment') == 'yes'
                and h.get('clip_cause') != 'window'),
            ','.join(fams) or '-',
            'yes' if czc else 'no', 'yes' if cytb else 'no',
            'yes' if 'hok_gef' in fams else 'no',
            off.get('czc', ''), off.get('hok_gef', ''),
            # recorded as context, NEVER as a membership test
            'yes' if (czc and cytb) else 'no',
            'no IS called' if not ish else '',
        ])

    with open(a.out, 'w') as fh:
        fh.write('\t'.join(hdr) + '\n')
        for r in rows:
            fh.write('\t'.join(str(v) for v in r) + '\n')

    # index by NAME -- these were positional and silently went stale when
    # columns were inserted ahead of them
    ix = {k: i for i, k in enumerate(hdr)}
    c = {}
    for r in rows:
        c[r[ix['vehicle_class']]] = c.get(r[ix['vehicle_class']], 0) + 1
    print(f'{len(rows)} loci -> {a.out}')
    for k in sorted(c):
        print(f'  {k:<11} {c[k]}')
    m = sum(1 for r in rows if r[ix['metal_module']] == 'yes')
    print(f'  with the czc+cytb metal module: {m}/{len(rows)}')
    nclip = sum(1 for r in rows if r[ix['vehicle_edge_clipped']] != '-')
    if nclip:
        print(f'  NOTE {nclip} locus/loci have a vehicle end cut by the window '
              f'edge -- coverage there is unknown, not decayed')


if __name__ == '__main__':
    main()
