#!/usr/bin/env python3
"""
Extract +/-20 kb windows around TADB3 relBE loci, located by LOCUS TAG.

Unlike the Stage 1-3 workflow (which finds relBE by blastn of the 527 bp
Enterobacteriaceae seed), TADB3 relBE spans distant phyla where that seed does
not hit. Loci are therefore located by the locus tags given in the TADB3 table.

Conventions enforced (see CLAUDE.md):
  * window is flipped so relB is on the + strand
  * +/-20 kb, clamped at replicon ends, ACTUAL flank lengths recorded
  * one record per copy: {accession}__{contig}__{ta_id}__{strand}
  * truncation is recorded, never a reason to drop a locus

usage: extract_tadb3_loci.py <parsed.tsv> <gbk_dir> <out_loci_dir> <out_table.tsv>
       [overrides.tsv]

Records lacking /locus_tag qualifiers are handled via an optional overrides
table giving explicit coordinates (see data/tadb3/locus_overrides.tsv).
"""
import os, re, sys

FLANK = int(os.environ.get('RELBE_FLANK', 15000))   # +/-15 kb: see PIPELINE.md F2


def parse_gbk(path):
    """Return (accession, length, topology, seq, [features])."""
    txt = open(path).read()
    loc = txt.split('\n', 1)[0]
    parts = loc.split()
    length = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else 0
    topology = 'circular' if 'circular' in loc else 'linear'
    acc = parts[1] if len(parts) > 1 else os.path.basename(path).split('.')[0]

    # sequence
    seq = ''
    if '\nORIGIN' in txt:
        origin = txt.split('\nORIGIN', 1)[1]
        seq = ''.join(re.findall(r'[acgtnACGTN]+', origin.split('//')[0].replace('\n', ' ')))
        seq = ''.join(c for c in seq if c.upper() in 'ACGTN').upper()

    # features
    feats = []
    if '\nFEATURES' not in txt:
        return acc, length, topology, seq, feats
    fb = txt.split('\nFEATURES', 1)[1].split('\nORIGIN', 1)[0]
    cur = None
    last_q = None
    for ln in fb.split('\n'):
        m = re.match(r'^ {5}(\S+)\s+(.+)$', ln)
        if m:
            if cur:
                feats.append(cur)
            cur = {'type': m.group(1), 'loc': m.group(2), 'q': {}}
            last_q = None
        elif cur is not None:
            s = ln.strip()
            if s.startswith('/'):
                mm = re.match(r'/(\w+)=?"?([^"]*)"?', s)
                if mm:
                    cur['q'][mm.group(1)] = mm.group(2)
                    last_q = mm.group(1)
            elif last_q and cur['loc'].count('(') > cur['loc'].count(')'):
                cur['loc'] += s          # continued location
            elif last_q:
                cur['q'][last_q] += ' ' + s
    if cur:
        feats.append(cur)
    return acc, length, topology, seq, feats


def find_tag(feats, tag):
    """Locate a feature by locus_tag / old_locus_tag / gene / protein_id."""
    for key in ('locus_tag', 'old_locus_tag', 'gene', 'protein_id'):
        for f in feats:
            if f['type'] not in ('CDS', 'gene'):
                continue
            v = f['q'].get(key, '')
            if v and tag.lower() in [x.strip().lower() for x in v.split(',')]:
                return f, key
    return None, None


def span(f):
    nums = [int(x) for x in re.findall(r'(\d+)', f['loc'])]
    return (min(nums), max(nums), 'complement' in f['loc']) if nums else None


def rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def load_overrides(path):
    ov = {}
    if not path or not os.path.exists(path):
        return ov
    for ln in open(path):
        if ln.startswith('#') or not ln.strip():
            continue
        f = ln.rstrip('\n').split('\t')
        ov[f[0]] = dict(acc=f[1], es=int(f[2]), ee=int(f[3]),
                        bs=int(f[4]), be=int(f[5]), bc=f[6].lower() == 'yes')
    return ov


def main():
    parsed, gbkdir, outdir, outtab = sys.argv[1:5]
    overrides = load_overrides(sys.argv[5] if len(sys.argv) > 5 else None)
    os.makedirs(outdir, exist_ok=True)
    rows = [l.rstrip('\n').split('\t') for l in open(parsed)]
    hdr, rows = rows[0], rows[1:]
    ci = {c: i for i, c in enumerate(hdr)}

    cache = {}
    out = []
    seen = set()
    for r in rows:
        acc = r[ci['accession']]
        ta = r[ci['ta_id']]
        key = (acc, r[ci['relE_locus_tag']], r[ci['relB_locus_tag']])
        if key in seen:
            print(f'  {ta:<9} {acc:<12} DUPLICATE locus, skipped')
            continue
        seen.add(key)

        path = os.path.join(gbkdir, acc + '.gbk')
        if not os.path.exists(path):
            print(f'  {ta:<9} {acc:<12} NO GBK FILE')
            out.append([ta, acc, '', '', '', '', '', '', '', '', '', 'missing_gbk'])
            continue
        if acc not in cache:
            cache[acc] = parse_gbk(path)
        gacc, glen, topo, seq, feats = cache[acc]

        if ta in overrides:
            o = overrides[ta]
            se, ee, ce = o['es'], o['ee'], False
            sb, eb, cb = o['bs'], o['be'], o['bc']
            fe = fb = True
        else:
            fe, ke = find_tag(feats, r[ci['relE_locus_tag']])
            fb, kb = find_tag(feats, r[ci['relB_locus_tag']])
        if not fe or not fb:
            miss = []
            if not fe: miss.append('relE:' + r[ci['relE_locus_tag']])
            if not fb: miss.append('relB:' + r[ci['relB_locus_tag']])
            print(f'  {ta:<9} {acc:<12} TAG NOT FOUND  {" ".join(miss)}')
            out.append([ta, acc, '', '', '', '', '', '', '', '', '', 'tag_not_found:' + ';'.join(miss)])
            continue

        if ta not in overrides:
            se, ee, ce = span(fe)
            sb, eb, cb = span(fb)
        lo, hi = min(se, sb), max(ee, eb)
        strand = 'minus' if cb else 'plus'     # orientation defined by relB

        ws, we = max(1, lo - FLANK), min(len(seq), hi + FLANK)
        up_actual = lo - ws
        dn_actual = we - hi
        trunc = []
        if up_actual < FLANK: trunc.append('left')
        if dn_actual < FLANK: trunc.append('right')

        win = seq[ws - 1:we]
        if strand == 'minus':
            win = rc(win)
            rel_s, rel_e = we - hi + 1, we - lo + 1
            up_final, dn_final = dn_actual, up_actual
        else:
            rel_s, rel_e = lo - ws + 1, hi - ws + 1
            up_final, dn_final = up_actual, dn_actual

        lid = f'{acc}__{gacc}__{ta}__{strand}'
        fp = os.path.join(outdir, lid + '.fa')
        with open(fp, 'w') as fh:
            fh.write(f'>{lid} orig={gacc}:{ws}-{we} strand={strand} '
                     f'flank_up={up_final} flank_dn={dn_final} relBE={rel_s}-{rel_e} '
                     f'truncated={",".join(trunc) or "no"}\n')
            for i in range(0, len(win), 70):
                fh.write(win[i:i + 70] + '\n')

        out.append([ta, acc, gacc, topo, str(glen), str(ws), str(we), str(len(win)),
                    str(up_final), str(dn_final), f'{rel_s}-{rel_e}',
                    ','.join(trunc) or 'no'])
        flag = '  TRUNCATED:' + ','.join(trunc) if trunc else ''
        print(f'  {ta:<9} {acc:<12} {strand:<5} relBE {hi-lo+1:>5} bp  '
              f'win {len(win):>6} bp  up {up_final:>6} dn {dn_final:>6}{flag}')

    with open(outtab, 'w') as fh:
        fh.write('ta_id\taccession\tcontig\ttopology\treplicon_len\twindow_start\t'
                 'window_end\twindow_len\tflank_up\tflank_dn\trelBE_in_window\ttruncated\n')
        for o in out:
            fh.write('\t'.join(o) + '\n')
    ok = sum(1 for o in out if o[2])
    print(f'\n{ok}/{len(out)} loci extracted -> {outdir}')


if __name__ == '__main__':
    main()
