#!/usr/bin/env python3
"""
F4 layer 2 + F5 -- IS elements by homology, with tiling and nesting.

  is_search.py --locus loci/plasmids/X.fa --db db/isfinder/ISfinder \
               --out results/is/X.tsv --bait-start 20001

NO LENGTH CUTOFF, anywhere. The whole Stage 5 test depends on a ~64 bp ISEsa1
remnant (cf. the ~75 bp ISAba125 remnant that identified Tn125 for blaNDM), so
--min-length defaults to 0 and is asserted to stay there.

Tiling uses TOLERANCE=30, never 0. PIPELINE.md records a real incident where a
3 bp overlap at zero tolerance discarded a 1,540 bp ISEsa1 alignment, making a
near-complete element look like a bare stub and inverting the biological
reading. Adjacent elements legitimately share a few bases at TSDs.

A low query_cov has two unrelated causes and the Stage 5 test only cares about
one of them. An element can be a genuine decayed remnant sitting in open DNA,
or it can be whole in the real replicon and merely sliced by the edge of the
+/-20 kb window we chose. Both come out of blastn as `is_fragment=yes`. Counting
them together would put phantom remnants at exactly +/-20000 in the ISEsa1
offset distribution -- the one plot the ancestral-vehicle conclusion rests on.
So edge contact is recorded explicitly, with the two causes kept apart:

  clipped_at_window_edge  no | left | right | both
  clip_cause              -  | window   (we cut it: discount as evidence)
                             | replicon (the replicon/contig truly ends here:
                                         real truncation, a signal per CLAUDE.md)
                             | mixed    (one of each)
  clipped_missing_bp      how much of the element lies beyond the edge

Nothing is filtered. `is_fragment` keeps its plain meaning (query_cov < 50) and
the new columns say whether that number can be read as decay.
"""

import argparse, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolpath import resolve  # noqa: E402

FMT = ('sseqid pident length qstart qend sstart send evalue bitscore slen')
TOLERANCE = 30


def hdr_field(path, key, default=None):
    with open(path) as fh:
        h = fh.readline()
    m = re.search(rf'{key}=(\S+)', h)
    return m.group(1) if m else default


def seq_len(path):
    with open(path) as fh:
        fh.readline()                      # header
        return sum(len(l.strip()) for l in fh)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--locus', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--min-length', type=int, default=0)
    p.add_argument('--evalue', default='1e-3')
    p.add_argument('--max-target', type=int, default=10000)
    p.add_argument('--tolerance', type=int, default=TOLERANCE)
    p.add_argument('--blast-bin',
                   default=os.environ.get('BLAST_BIN', ''),
                   help='dir holding blastn; empty => resolve from PATH')
    a = p.parse_args()

    assert a.min_length == 0, (
        'min_length must stay 0: the ~64 bp ISEsa1 remnant is the key evidence '
        'for the ancestral-vehicle test (PIPELINE.md section 3)')
    assert a.tolerance > 0, 'tolerance must be > 0 (PIPELINE.md section 4 trap)'

    lid = os.path.basename(a.locus).rsplit('.', 1)[0]
    bait = hdr_field(a.locus, 'bait', '1-1')
    bait_start = int(bait.split('-')[0])
    win_len = seq_len(a.locus)
    # post-reorientation, so truncated_up is the LEFT edge of this FASTA
    trunc_left = hdr_field(a.locus, 'truncated_up', 'no') == 'yes'
    trunc_right = hdr_field(a.locus, 'truncated_dn', 'no') == 'yes'

    cmd = [resolve('blastn', a.blast_bin),
           '-query', a.locus, '-db', a.db, '-task', 'blastn',
           '-evalue', a.evalue, '-outfmt', '6 ' + FMT,
           '-max_target_seqs', str(a.max_target),
           '-dust', 'no', '-soft_masking', 'false']
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit('blastn failed:\n' + r.stderr)

    hsps = []
    for line in r.stdout.splitlines():
        f = line.split('\t')
        d = dict(zip(FMT.split(), f))
        qs, qe = sorted((int(d['qstart']), int(d['qend'])))
        ss, se = int(d['sstart']), int(d['send'])
        name = d['sseqid']
        parts = name.split('_')
        hsps.append({
            'element': parts[0],
            'family': parts[1] if len(parts) > 1 else '',
            'group': parts[2] if len(parts) > 2 else '',
            'sseqid': name,
            'pident': float(d['pident']),
            'length': int(d['length']),
            'win_start': qs, 'win_end': qe,
            'sub_start': min(ss, se), 'sub_end': max(ss, se),
            'strand': '-' if ss > se else '+',
            'slen': int(d['slen']),
            'bitscore': float(d['bitscore']),
            'evalue': d['evalue'],
        })

    # ---- tiling: keep best-scoring HSPs, allow TOLERANCE bp of overlap ----
    hsps.sort(key=lambda h: -h['bitscore'])
    kept = []
    for h in hsps:
        clash = False
        for k in kept:
            ov = min(h['win_end'], k['win_end']) - max(h['win_start'], k['win_start']) + 1
            if ov > a.tolerance:
                clash = True
                break
        if not clash:
            kept.append(h)
    kept.sort(key=lambda h: h['win_start'])

    # ---- coverage class; fragments are KEPT, never filtered ----
    for h in kept:
        cov = 100.0 * h['length'] / h['slen']
        h['query_cov'] = round(cov, 1)
        h['cov_class'] = ('complete' if cov >= 90 else
                          'partial' if cov >= 50 else 'fragment')
        h['is_fragment'] = 'yes' if cov < 50 else 'no'
        h['offset_from_bait'] = h['win_start'] - bait_start

    # ---- window-edge clipping ----
    # An alignment touching an edge is only CLIPPED if part of the subject is
    # left unaccounted for on that same side -- i.e. the element carries on past
    # where we stopped looking. Which subject end that is depends on strand:
    # on '-' the subject runs backwards relative to the window.
    for h in kept:
        if h['strand'] == '+':
            miss_l, miss_r = h['sub_start'] - 1, h['slen'] - h['sub_end']
        else:
            miss_l, miss_r = h['slen'] - h['sub_end'], h['sub_start'] - 1
        sides, causes, missing = [], [], 0
        if h['win_start'] <= 1 and miss_l > 0:
            sides.append('left')
            causes.append('replicon' if trunc_left else 'window')
            missing += miss_l
        if h['win_end'] >= win_len and miss_r > 0:
            sides.append('right')
            causes.append('replicon' if trunc_right else 'window')
            missing += miss_r
        h['clipped_at_window_edge'] = ('both' if len(sides) == 2
                                       else sides[0] if sides else 'no')
        h['clip_cause'] = ('-' if not causes else
                           causes[0] if len(set(causes)) == 1 else 'mixed')
        h['clipped_missing_bp'] = missing

    # ---- nesting: element B lies inside the window span of element A, where A
    # survives as two HSPs whose SUBJECT ranges are near-contiguous ----
    by_elem = {}
    for h in kept:
        by_elem.setdefault(h['element'], []).append(h)
    for h in kept:
        h['nested_in'] = ''
        h['interrupts'] = ''
    for elem, parts in by_elem.items():
        if len(parts) < 2:
            continue
        # PIPELINE.md section 4: element A counts as SPLIT -- and therefore as
        # something another element can be nested inside -- only if its HSPs are
        # near-contiguous in SUBJECT coordinates, i.e. they are two pieces of ONE
        # copy (ISEsa1 bases 1-1540 and 1533-1596). Without this check, several
        # independent FULL-LENGTH copies of the same element are misread as one
        # split copy, and everything between the outermost two is falsely called
        # nested. That fires hard on IS26/IS15DI arrays, which are near-identical
        # and recur many times per plasmid.
        bysub = sorted(parts, key=lambda x: x['sub_start'])
        contiguous = all(
            abs(bysub[i + 1]['sub_start'] - bysub[i]['sub_end']) <= a.tolerance
            for i in range(len(bysub) - 1))
        total_sub = sum(x['sub_end'] - x['sub_start'] + 1 for x in bysub)
        one_copy = total_sub <= 1.1 * bysub[0]['slen']
        if not (contiguous and one_copy):
            continue
        parts.sort(key=lambda x: x['win_start'])
        lo, hi = parts[0]['win_end'], parts[-1]['win_start']
        for h in kept:
            if h['element'] == elem:
                continue
            if h['win_start'] >= lo - a.tolerance and h['win_end'] <= hi + a.tolerance:
                h['nested_in'] = elem
                h['interrupts'] = elem

    cols = ['locus_id', 'element', 'family', 'group', 'win_start', 'win_end',
            'length', 'strand', 'pident', 'query_cov', 'cov_class',
            'is_fragment', 'clipped_at_window_edge', 'clip_cause',
            'clipped_missing_bp', 'offset_from_bait', 'sub_start', 'sub_end',
            'slen', 'nested_in', 'interrupts', 'evalue', 'bitscore']
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for h in kept:
            h['locus_id'] = lid
            fh.write('\t'.join(str(h[c]) for c in cols) + '\n')

    nf = sum(1 for h in kept if h['is_fragment'] == 'yes')
    nclip = sum(1 for h in kept
                if h['is_fragment'] == 'yes' and h['clip_cause'] == 'window')
    print(f'{lid}: {len(hsps)} HSPs -> {len(kept)} tiled '
          f'({nf} fragment, of which {nclip} cut by the window edge; '
          f'{nf - nclip} usable as remnants) -> {a.out}')


if __name__ == '__main__':
    main()
