#!/usr/bin/env python3
"""
F2 -- cut the +/-flank window around each selected locus and reorient.

  extract_loci.py --selected results/loci_selected.tsv --fasta-dir data/public_16/fasta \
                  --out-dir loci --out-table results/loci.tsv --flank 20000

Locating is already done (Tier A/B); this step only cuts. The actual cutting
lives in scripts/relbe_window.py, which handles circular wrap-around, the
whole-replicon case, and keeps truncation flags in the same coordinate frame as
the flanks -- all three of which the previous extractor got wrong.

Reorientation matters more than it looks: relBE sits on the minus strand in
some replicons and the plus strand in others. Without a common orientation
every downstream offset is mirrored for half the loci and nothing lines up.
"""

import argparse, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from relbe_window import extract_window          # noqa: E402
from gbk_tools import read_fasta, write_tsv      # noqa: E402
from toolpath import resolve                     # noqa: E402

HDR = ['locus_id', 'accession', 'contig', 'organism', 'genus', 'replicon_type',
       'replicon_name', 'replicon_len', 'topology', 'is_complete',
       'hit_start_orig', 'hit_end_orig', 'strand_orig', 'reoriented',
       'window_len', 'flank_up', 'flank_dn', 'bait_in_window',
       'truncated_up', 'truncated_dn', 'wrapped', 'whole_replicon',
       'pident', 'qcov', 'source_interval']


def seq_from_blastdb(acc, db, blast_bin=''):
    """Pull a full record out of a BLAST database by accession.

    This is what makes a large sweep tractable. A pre-formatted database from
    `update_blastdb.pl` (nt_prok and friends) ships no FASTAs, but it already
    CONTAINS every sequence -- so the flanking DNA for a hit is already local
    and needs no second download. Requires the db to carry seq-ids, which is
    true of NCBI's pre-formatted sets and of anything built with -parse_seqids.

    The whole record is fetched rather than a range, because a circular window
    may wrap the origin and the wrap arithmetic lives in relbe_window.
    """
    for entry in (acc, acc.split('.')[0]):
        r = subprocess.run([resolve('blastdbcmd', blast_bin), '-db', db,
                            '-entry', entry, '-outfmt', '%s'],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip().replace('\n', '').upper()
    return None


def find_fasta(d, acc):
    for cand in (acc, acc.split('.')[0]):
        for ext in ('.fna', '.fa', '.fasta'):
            p = os.path.join(d, cand + ext)
            if os.path.exists(p):
                return p
        import glob
        g = glob.glob(os.path.join(d, cand + '.*' + ext.lstrip('.')))
        if g:
            return g[0]
    import glob
    g = glob.glob(os.path.join(d, acc.split('.')[0] + '*'))
    g = [x for x in g if x.endswith(('.fna', '.fa', '.fasta'))]
    return g[0] if g else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--selected', required=True)
    p.add_argument('--fasta-dir', default='',
                   help='directory of source FASTAs (curated sets)')
    p.add_argument('--blastdb', default='',
                   help='BLAST db to pull sequence from instead; required for '
                        'pre-formatted sets like nt_prok, which ship no FASTAs')
    p.add_argument('--blast-bin', default=os.environ.get('BLAST_BIN', ''))
    p.add_argument('--out-dir', required=True)
    p.add_argument('--out-table', required=True)
    p.add_argument('--flank', type=int, default=20000)
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    with open(a.selected) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        rows = [dict(zip(hdr, l.rstrip('\n').split('\t')))
                for l in fh if l.strip()]
    rows = [r for r in rows if r.get('selected') == 'yes']

    out, nwrap, ntrunc = [], 0, 0
    if not a.fasta_dir and not a.blastdb:
        sys.exit('give --fasta-dir or --blastdb')

    nmiss = 0
    for r in rows:
        seq = None
        if a.fasta_dir:
            fp = find_fasta(a.fasta_dir, r['accession'])
            if fp:
                _, seq = read_fasta(fp)
        if seq is None and a.blastdb:
            seq = seq_from_blastdb(r['accession'], a.blastdb, a.blast_bin)
        if seq is None:
            nmiss += 1
            print(f"  SKIP {r['accession']}: not in --fasta-dir or --blastdb")
            continue
        w = extract_window(seq, int(r['hit_start']), int(r['hit_end']),
                           r['strand'] == 'minus', a.flank,
                           r['topology'], len(seq))

        lid = r['locus_id']
        with open(os.path.join(a.out_dir, lid + '.fa'), 'w') as fo:
            fo.write(f">{lid} orig={r['accession']}:{w['orig']} "
                     f"strand={w['strand']} flank_up={w['flank_up']} "
                     f"flank_dn={w['flank_dn']} "
                     f"bait={w['rel_start']}-{w['rel_end']} "
                     f"truncated_up={'yes' if w['truncated_up'] else 'no'} "
                     f"truncated_dn={'yes' if w['truncated_dn'] else 'no'} "
                     f"wrapped={'yes' if w['wrapped'] else 'no'} "
                     f"whole_replicon={'yes' if w['whole_replicon'] else 'no'}\n")
            s = w['seq_out']
            for i in range(0, len(s), 70):
                fo.write(s[i:i + 70] + '\n')

        nwrap += bool(w['wrapped'])
        ntrunc += bool(w['truncated_up'] or w['truncated_dn'])
        out.append([lid, r['accession'], r['accession'], r['organism'],
                    r['genus'], r['replicon_type'], r['replicon_name'],
                    r['replicon_len'], r['topology'], r['is_complete'],
                    r['hit_start'], r['hit_end'], r['strand'],
                    'revcomp' if w['strand'] == 'minus' else 'none',
                    w['window_len'], w['flank_up'], w['flank_dn'],
                    f"{w['rel_start']}-{w['rel_end']}",
                    'yes' if w['truncated_up'] else 'no',
                    'yes' if w['truncated_dn'] else 'no',
                    'yes' if w['wrapped'] else 'no',
                    'yes' if w['whole_replicon'] else 'no',
                    r['pident'], r['qcov'], w['orig']])

    write_tsv(a.out_table, HDR, out)
    print(f'{len(out)} windows -> {a.out_dir}/  table -> {a.out_table}')
    print(f'  wrapped across origin: {nwrap}   short flank (kept): {ntrunc}')
    if nmiss:
        print(f'  WARNING {nmiss} selected locus/loci had no retrievable sequence')


if __name__ == '__main__':
    main()
