#!/usr/bin/env python3
"""
Stage 3/5 -- find the direct repeat (target-site duplication) flanking each
IS call, and compare it across loci: CLAUDE.md's "distinct DRs across loci =
independent events".

  relbe_junctions.py --is run_all/is --loci run_all/loci \
                     --out results/junctions.tsv

HOW A TSD IS FOUND
-------------------
A classic DDE transposase makes a staggered double-strand cut in the target
DNA; when the element is ligated in and the gaps are filled, the few bp of
host sequence that were in the staggered overhang end up duplicated -- once
immediately upstream of the insertion, once immediately downstream. A TSD is
therefore ONE comparison per element, not two independent ones: take the N bp
of window sequence immediately before the hit and the N bp immediately after
it, and look for the longest exact match between a SUFFIX of the upstream
flank and a PREFIX of the downstream flank.

A v1 of this script got exactly this wrong -- it looped over 'left side' and
'right side' as if they were separate searches, but a TSD search only makes
sense as a single two-sided comparison, so both loop iterations ended up
comparing the SAME pair of flanks and printing the identical result twice
under two different labels. That looked like independent confirmation from
two boundaries; it was one comparison counted twice. Fixed here: one row per
IS hit, one comparison, using BOTH flanks together.

A BOUNDARY IS ONLY SEARCHED IF IT IS GENUINELY FLANKED BY NEUTRAL SEQUENCE
------------------------------------------------------------------------------
Three ways a boundary can fail to have real host flank immediately next to
it, each checked and reported rather than assumed:

  - the window itself cuts it off (`clipped_at_window_edge`) -- there is no
    sequence on that side within the window to compare;
  - it touches another IS call (within `--flank` bp) -- what is "outside"
    this element is in fact a DIFFERENT element, not host DNA. This is the
    ISEsa1/ISKpn21 case in CP018961: ISEsa1 is split into a 1540 bp piece and
    a 64 bp remnant by ISKpn21 landing inside it, and the boundary between
    them is an interruption, not an insertion terminus;
  - it is a fragment, and the window boundary is not where the reference
    sequence itself ends (`sub_start`/`sub_end` vs `slen`, strand-aware) --
    i.e. this is simply where the BLAST alignment stopped covering a
    degraded copy, not the original element's true edge. What lies beyond it
    could be more unrecognised IS sequence, not host DNA.

Both flanks have to be clean for the comparison to mean anything, so a hit
with even one bad boundary is excluded outright, with the reason recorded.

NOT EVERY IS FAMILY MAKES A TSD AT ALL -- AND THIS MATTERS FOR OUR VEHICLE
------------------------------------------------------------------------------
IS110 (serine-recombinase-like, site-specific) and IS481 do not generate TSDs
by the staggered-cut mechanism above; a conserved sequence flanking an IS110
hit is more likely a conserved target motif the recombinase recognises (or
simply conserved host sequence the element happens to sit in) than a
duplication signature. ISEsa2 IS IS110. A clean, repeatable TSD is therefore
only informative for ISEsa1 (ISNCY); for ISEsa2 the script still reports
whatever it finds (`tsd_expected=no`) but that finding cannot distinguish one
ancestral insertion from several independent ones landing on the same
targeted site, and should not be read as if it could.

CROSS-LOCUS COMPARISON IS THE POINT
---------------------------------------
For each element name, loci sharing an EXACT TSD are grouped (`same_dr_group`,
content-addressed like `structural_groups.py`, so colour is stable). Loci
whose TSD differs, or where none was found, are left ungrouped.
"""

import argparse, csv, glob, hashlib, os, sys
from collections import defaultdict

SPARE = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02',
         '#a6761d', '#666666', '#1f78b4', '#b2182b']

TSD_BY_FAMILY = {'IS110': 'no', 'IS481': 'no', 'IS200/IS605': 'no'}


def dr_colour(seq):
    h = int(hashlib.sha1(seq.encode()).hexdigest(), 16)
    return SPARE[h % len(SPARE)]


def read_fasta_one(path):
    seq = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            if line.startswith('>'):
                break
            seq.append(line.strip())
    return ''.join(seq)


def best_match(left, right, min_k, max_k, min_id):
    profile = []
    best_k, best_id, best_seq = 0, 0.0, ''
    for k in range(min_k, min(max_k, len(left), len(right)) + 1):
        a, b = left[-k:], right[:k]
        pid = sum(1 for x, y in zip(a, b) if x.upper() == y.upper()) / k
        profile.append((k, round(pid * 100, 1)))
        if pid >= min_id:
            best_k, best_id, best_seq = k, pid, b
    return best_k, best_id, best_seq, profile


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--is', dest='is_dir', required=True)
    p.add_argument('--loci', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--out-profile', default='')
    p.add_argument('--flank', type=int, default=25,
                   help='bp of flank searched on each side, per CLAUDE.md')
    p.add_argument('--min-k', type=int, default=2)
    p.add_argument('--max-k', type=int, default=20)
    p.add_argument('--min-id', type=float, default=0.90)
    p.add_argument('--elements', default='',
                   help='restrict to these element names (comma-separated); '
                        'empty = every IS call')
    a = p.parse_args()

    files = sorted(glob.glob(os.path.join(a.is_dir, '*.is.tsv')))
    if not files:
        sys.exit(f'no *.is.tsv under {a.is_dir}')
    want = {e.strip() for e in a.elements.split(',') if e.strip()}

    rows_out, profile_out = [], []
    for f in files:
        lid = os.path.basename(f)[:-len('.is.tsv')]
        fa = os.path.join(a.loci, lid + '.fa')
        if not os.path.exists(fa):
            sys.exit(f'no window FASTA for {lid}: expected {fa}')
        win = read_fasta_one(fa)
        hits = list(csv.DictReader(open(f), delimiter='\t'))

        def touches_another(ws, we, side):
            # Does ANOTHER is call sit within --flank bp of this boundary?
            # If so what is "outside" this element here is a different
            # element, not neutral host sequence.
            for o in hits:
                if o is h:
                    continue
                ows, owe = int(o['win_start']), int(o['win_end'])
                if side == 'left' and 0 <= ws - owe <= a.flank:
                    return o['element']
                if side == 'right' and 0 <= ows - we <= a.flank:
                    return o['element']
            return ''

        for h in hits:
            if want and h['element'] not in want:
                continue
            ws, we = int(h['win_start']), int(h['win_end'])
            strand = h['strand']
            sub_start, sub_end, slen = (int(h['sub_start']), int(h['sub_end']),
                                        int(h['slen']))
            clip = h.get('clipped_at_window_edge', 'no')
            frag = h.get('cov_class') == 'fragment'

            # a few bp of slack: a split element's two HSPs can overlap or
            # gap slightly at the junction (TOLERANCE=30 elsewhere in this
            # pipeline, is_search.py), so exact equality under-recognises a
            # piece that genuinely does reach the reference's true end
            TOL = 5
            left_is_outer = (sub_start <= 1 + TOL) if strand == '+' else (sub_end >= slen - TOL)
            right_is_outer = (sub_end >= slen - TOL) if strand == '+' else (sub_start <= 1 + TOL)

            excl = []
            if clip in ('left', 'both'):
                excl.append('window_edge:left')
            if clip in ('right', 'both'):
                excl.append('window_edge:right')
            touch_l = touches_another(ws, we, 'left')
            if touch_l:
                excl.append(f'touches:{touch_l}:left')
            touch_r = touches_another(ws, we, 'right')
            if touch_r:
                excl.append(f'touches:{touch_r}:right')
            if frag and not left_is_outer:
                excl.append('fragment_internal:left')
            if frag and not right_is_outer:
                excl.append('fragment_internal:right')

            row = {'locus_id': lid, 'element': h['element'],
                  'family': h['family'], 'win_start': ws, 'win_end': we,
                  'strand': strand, 'cov_class': h['cov_class'],
                  'tsd_expected': TSD_BY_FAMILY.get(h['family'], 'yes')}

            if excl:
                row.update(excluded_reason=';'.join(excl), dr_seq='', dr_len=0,
                          dr_id='', left_flank='', right_flank='')
                rows_out.append(row)
                continue

            lo = max(0, ws - 1 - a.flank)
            hi = min(len(win), we + a.flank)
            left_flank = win[lo:ws - 1]
            right_flank = win[we:hi]
            k, pid, seq, profile = best_match(left_flank, right_flank,
                                              a.min_k, a.max_k, a.min_id)
            row.update(excluded_reason='', dr_seq=seq, dr_len=k,
                      dr_id=f'{pid * 100:.0f}' if k else '',
                      left_flank=left_flank, right_flank=right_flank)
            rows_out.append(row)
            for kk, idp in profile:
                profile_out.append({'locus_id': lid, 'element': h['element'],
                                    'k': kk, 'pct_id': idp})

    by_elem_dr = defaultdict(set)
    for r in rows_out:
        if r['dr_seq']:
            by_elem_dr[(r['element'], r['dr_seq'].upper())].add(r['locus_id'])
    group_of = {}
    for (elem, seq), loci in by_elem_dr.items():
        if len(loci) > 1:
            gid = 'DR' + hashlib.sha1(f'{elem}|{seq.upper()}'.encode()).hexdigest()[:6]
            for l in loci:
                group_of[(l, elem)] = (gid, dr_colour(seq.upper()), len(loci))

    cols = ['locus_id', 'element', 'family', 'win_start', 'win_end', 'strand',
            'cov_class', 'tsd_expected', 'dr_seq', 'dr_len', 'dr_id',
            'same_dr_group', 'same_dr_group_size', 'same_dr_colour',
            'left_flank', 'right_flank', 'excluded_reason']
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows_out:
            g = group_of.get((r['locus_id'], r['element']))
            r['same_dr_group'] = g[0] if g else ''
            r['same_dr_colour'] = g[1] if g else ''
            r['same_dr_group_size'] = g[2] if g else ''
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

    if a.out_profile:
        with open(a.out_profile, 'w') as fh:
            fh.write('locus_id\telement\tk\tpct_id\n')
            for r in profile_out:
                fh.write(f"{r['locus_id']}\t{r['element']}\t{r['k']}\t{r['pct_id']}\n")

    searched = [r for r in rows_out if not r['excluded_reason']]
    found = [r for r in searched if r['dr_seq']]
    print(f'{len(rows_out)} IS calls examined')
    print(f'  {len(rows_out) - len(searched)} excluded (window edge / '
          f'touches another IS / fragment-internal boundary)')
    print(f'  {len(searched)} searched, {len(found)} with a candidate TSD '
          f'(>= {a.min_id:.0%} id, {a.min_k}-{a.max_k} bp)')

    print(f'\ncross-locus agreement, by element:')
    by_elem = defaultdict(list)
    for r in found:
        by_elem[r['element']].append(r)
    for elem, rs in sorted(by_elem.items()):
        seqs = defaultdict(list)
        for r in rs:
            seqs[r['dr_seq'].upper()].append(r['locus_id'].split('__')[0])
        exp = TSD_BY_FAMILY.get(rs[0]['family'], 'yes')
        warn = ('  (not mechanistically expected -- site-specific family, '
                'see docstring)' if exp == 'no' else '')
        print(f'  {elem}{warn}')
        for seq, loci in sorted(seqs.items(), key=lambda kv: -len(kv[1])):
            tag = '  ** shared **' if len(loci) > 1 else ''
            print(f'    {seq:<16} {len(loci)} loci: {", ".join(loci)}{tag}')

    print(f'\n-> {a.out}')
    if a.out_profile:
        print(f'-> {a.out_profile}')


if __name__ == '__main__':
    main()
