#!/usr/bin/env python3
"""
One readable account of a TAflank run: what went in, what each filter removed,
and what came out at every layer.

  run_report.py --results results --out results/run_report.txt \
                --tsv results/run_report.tsv

Reads only the published tables, so the report is reproducible from results
alone and does not need the pipeline to be re-run. Anything absent is reported
as such rather than silently omitted -- a layer that did not run should be
visible, not invisible.
"""

import argparse, datetime, glob, os


def read_tsv(p):
    if not p or not os.path.exists(p):
        return None
    with open(p) as fh:
        lines = [l for l in fh if not l.startswith('#') and l.strip()]
    if not lines:
        return []
    hdr = lines[0].rstrip('\n').split('\t')
    return [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in lines[1:]]


def count_blast(p):
    if not os.path.exists(p):
        return None
    return sum(1 for l in open(p) if l.strip() and not l.startswith('#'))


def tally(rows, key):
    c = {}
    for r in rows or []:
        c[r.get(key, '')] = c.get(r.get(key, ''), 0) + 1
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', default='results')
    ap.add_argument('--out', default='')
    ap.add_argument('--tsv', default='')
    a = ap.parse_args()
    R = a.results
    J = lambda *p: os.path.join(R, *p)

    raw    = count_blast(J('blast', '01_raw.tsv'))
    tierA  = count_blast(J('blast', '02_tier_a.tsv'))
    tierB  = count_blast(J('blast', '03_selected.tsv'))
    hits   = read_tsv(J('bait_hits.tsv'))
    sel    = read_tsv(J('loci_selected.tsv'))
    loci   = read_tsv(J('loci.tsv'))
    mge    = read_tsv(J('mge_hits.tsv'))
    rep    = read_tsv(J('replicons.tsv'))
    ctx    = read_tsv(J('context.tsv'))
    audit  = read_tsv(J('filter_summary.tsv'))
    cargo  = []
    for f in glob.glob(J('cargo', '*.cargo.tsv')):
        cargo += read_tsv(f) or []
    blocks = []
    for f in glob.glob(J('blocks', '*.blocks.tsv')):
        blocks += read_tsv(f) or []

    L, M = [], []                      # text lines, (metric, value) pairs
    def add(label, value, note=''):
        L.append('  %-38s %10s  %s' % (label, value, note))
        M.append((label.strip(), value, note))
    def hdr(t):
        L.append(''); L.append(t); L.append('-' * len(t))

    params = ''
    for r in audit or []:
        if r.get('stage') == 'params':
            params = r.get('category', '')

    L.append('TAflank run report')
    L.append('=' * 60)
    L.append('  generated : %s' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
    L.append('  results   : %s' % os.path.abspath(R))
    if params:
        L.append('  parameters: %s' % params)

    hdr('DISCOVERY  (how the hits were narrowed)')
    if raw is None:
        add('blastn HSPs', 'n/a', 'results/blast/ not written')
    else:
        add('blastn HSPs, unfiltered', raw)
        drop_a = (raw - tierA) if tierA is not None else None
        add('after Tier A (sequence)', tierA if tierA is not None else 'n/a',
            '-%d removed' % drop_a if drop_a else 'none removed')
        drop_b = (tierA - tierB) if (tierA is not None and tierB is not None) else None
        add('after Tier B (record)', tierB if tierB is not None else 'n/a',
            '-%d removed' % drop_b if drop_b else 'none removed')
    for k, v in sorted(tally(hits, 'tier_a').items()):
        if k:
            add('  Tier A class: %s' % k, v)
    for r in audit or []:
        if r.get('stage') == 'deselect_reason':
            add('  deselected: %s' % r['category'][:30], r['n'])

    hdr('WINDOWS')
    if loci is None:
        add('windows extracted', 'n/a')
    else:
        add('windows extracted', len(loci))
        tot = sum(int(r.get('window_len') or 0) for r in loci)
        add('total sequence', '%s bp' % format(tot, ','))
        full = sum(1 for r in loci
                   if r.get('truncated_up') == 'no' and r.get('truncated_dn') == 'no')
        add('full flanks both sides', '%d/%d' % (full, len(loci)))
        add('truncated (kept, recorded)', len(loci) - full)
        add('wrapped the origin', sum(1 for r in loci if r.get('wrapped') == 'yes'))
        for k, v in sorted(tally(loci, 'replicon_type').items()):
            add('  replicon_type: %s' % k, v)

    hdr('MGE LAYER  (ISfinder, no length cutoff)')
    if mge is None:
        add('IS calls', 'n/a')
    else:
        add('IS calls', len(mge))
        for k in ('complete', 'partial', 'fragment'):
            add('  %s' % k, sum(1 for r in mge if r.get('cov_class') == k))
        nclip = sum(1 for r in mge if r.get('clip_cause') == 'window')
        add('  of which cut by window edge', nclip)
        add('  fragments usable as remnants',
            sum(1 for r in mge if r.get('is_fragment') == 'yes'
                and r.get('clip_cause') != 'window'))
        add('distinct elements', len({r['element'] for r in mge}))
        add('nesting events', sum(1 for r in mge if r.get('nested_in')))
        add('loci with ISEsa1', len({r['locus_id'] for r in mge
                                     if r['element'] == 'ISEsa1'}))
        add('loci with ISEsa2', len({r['locus_id'] for r in mge
                                     if r['element'] == 'ISEsa2'}))

    hdr('CARGO')
    if not cargo:
        add('cargo hits', 'n/a' if loci is None else 0)
    else:
        add('cargo hits', len(cargo))
        n = len(loci or [])
        for fam in ('czc', 'cytb', 'hok_gef', 'tra', 'qnr', 'intI'):
            k = len({r['locus_id'] for r in cargo if r['family'] == fam})
            if k:
                add('  loci with %s' % fam, '%d/%d' % (k, n))

    hdr('CONTEXT  (vehicle class)')
    if ctx is None:
        add('classified', 'n/a')
    else:
        for k, v in sorted(tally(ctx, 'vehicle_class').items()):
            add(k, v)
        add('with czc+cytb metal module',
            sum(1 for r in ctx if r.get('metal_module') == 'yes'))

    hdr('REPLICON TYPING  (PlasmidFinder)')
    if rep is None:
        add('replicon hits', 'n/a')
    else:
        add('replicon hits', len(rep))
        add('records typed', len({r['accession'] for r in rep}))
        inc = tally(rep, 'replicon')
        for k, v in sorted(inc.items(), key=lambda x: -x[1])[:6]:
            add('  %s' % k, v)

    hdr('FEATURE TRACK  (figure input)')
    if not blocks:
        add('drawable features', 'n/a')
    else:
        add('drawable features', len(blocks))
        for k, v in sorted(tally(blocks, 'colour_source').items()):
            add('  coloured by %s' % k, v)
        srcs = tally(blocks, 'annot_source')
        for k, v in sorted(srcs.items()):
            add('  annotation source: %s' % (k or 'none'), v)
        empty = [l for l in (loci or [])
                 if not any(b['locus_id'] == l['locus_id'] and
                            b['feature_type'] == 'CDS' for b in blocks)]
        if empty:
            add('loci with NO gene track', len(empty),
                '<- needs Bakta: ' + ','.join(e['accession'] for e in empty))

    L.append('')
    text = '\n'.join(L)
    print(text)
    if a.out:
        os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
        open(a.out, 'w').write(text + '\n')
    if a.tsv:
        with open(a.tsv, 'w') as fh:
            fh.write('metric\tvalue\tnote\n')
            for k, v, n in M:
                fh.write('%s\t%s\t%s\n' % (k, v, n))


if __name__ == '__main__':
    main()
