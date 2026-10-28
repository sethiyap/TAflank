#!/usr/bin/env python3
"""
F4 layer 5 -- replicon (Inc) typing by PlasmidFinder, via plain blastn.

  plasmidfinder.py --replicon data/public_16/fasta/CP018961.1.fna \
                   --db db/plasmidfinder/plasmidfinder.fa --out results/pf/CP018961.tsv

No CGE tooling, no conda: the PlasmidFinder database is plain FASTA and the
native BLAST already on disk types it directly. Thresholds are Acman et al.
2022 (Nat Commun 13:1131): coverage >= 80%, identity >= 95%, and overlapping
hits on one replicon resolved in favour of the higher identity.

Run on the WHOLE PARENT REPLICON, never the +/-20 kb window: Inc rep genes are
rarely within 20 kb of relBE, so a window-only run returns empty for nearly
everything.

The replicon gene is the QUERY and the replicon the SUBJECT, so coverage is
qcovs -- coverage of the rep gene. Getting this backwards is the same
qlen-vs-slen trap that once misclassified all nine loci.
"""

import argparse, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolpath import resolve  # noqa: E402

FMT = 'qseqid qlen length pident sstart send evalue bitscore qcovs'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--replicon', default='',
                   help='FASTA of the whole parent replicon')
    p.add_argument('--blastdb', default='',
                   help='pull the replicon from this BLAST db instead; needed '
                        'for pre-formatted sets like nt_prok that ship no FASTAs')
    p.add_argument('--efetch', action='store_true',
                   help='fetch the whole replicon from NCBI; for hits imported '
                        'from a web BLAST, where no local copy exists')
    p.add_argument('--efetch-cache', default='data/efetch_cache')
    p.add_argument('--db', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--accession', default='')
    p.add_argument('--min-cov', type=float, default=80.0)
    p.add_argument('--min-id', type=float, default=95.0)
    p.add_argument('--blast-bin',
                   default=os.environ.get('BLAST_BIN', ''),
                   help='dir holding blastn; empty => resolve from PATH')
    a = p.parse_args()

    acc = a.accession or os.path.basename(a.replicon).split('.')[0]

    subject = a.replicon
    tmp = ''
    if not subject and a.blastdb:
        # Materialise the replicon from the db just long enough to BLAST it.
        if not acc:
            sys.exit('--blastdb needs --accession')
        r = subprocess.run([resolve('blastdbcmd', a.blast_bin), '-db', a.blastdb,
                            '-entry', acc, '-outfmt', '%f'],
                           capture_output=True, text=True)
        if r.returncode or not r.stdout.strip():
            r = subprocess.run([resolve('blastdbcmd', a.blast_bin), '-db', a.blastdb,
                                '-entry', acc.split('.')[0], '-outfmt', '%f'],
                               capture_output=True, text=True)
        if r.returncode or not r.stdout.strip():
            sys.exit(f'could not retrieve {acc} from {a.blastdb}')
        tmp = f'_{acc}.replicon.fa'
        with open(tmp, 'w') as fh:
            fh.write(r.stdout)
        subject = tmp
    if not subject and a.efetch:
        # The WHOLE record, not a window: Inc rep genes are rarely within
        # 20 kb of the bait, so a window-only search returns nothing useful.
        import urllib.parse, urllib.request, time
        os.makedirs(a.efetch_cache, exist_ok=True)
        cf = os.path.join(a.efetch_cache, f'{acc}_full.fa')
        if not os.path.exists(cf):
            q = {'db': 'nuccore', 'id': acc, 'rettype': 'fasta',
                 'retmode': 'text'}
            kf = os.path.expanduser('~/.ncbi_api_key')
            if os.path.exists(kf):
                q['api_key'] = open(kf).read().strip()
            url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
                   + urllib.parse.urlencode(q))
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url, timeout=180) as r:
                        txt = r.read().decode()
                    if txt.startswith('>'):
                        open(cf, 'w').write(txt)
                        break
                except Exception:
                    pass
                time.sleep(2 ** attempt)
        if os.path.exists(cf):
            subject = cf
    if not subject:
        sys.exit('give --replicon, --blastdb or --efetch')

    cmd = [resolve('blastn', a.blast_bin),
           '-query', a.db, '-subject', subject, '-task', 'blastn',
           '-evalue', '1e-10', '-outfmt', '6 ' + FMT,
           '-dust', 'no', '-soft_masking', 'false']
    r = subprocess.run(cmd, capture_output=True, text=True)
    if tmp and os.path.exists(tmp):
        os.remove(tmp)
    if r.returncode:
        sys.exit('blastn failed:\n' + r.stderr)

    hits = []
    for line in r.stdout.splitlines():
        f = line.split('\t')
        d = dict(zip(FMT.split(), f))
        cov, pid = float(d['qcovs']), float(d['pident'])
        if cov < a.min_cov or pid < a.min_id:
            continue
        q = d['qseqid']                       # "<group>|<replicon>_<n>_<acc>"
        grp, rest = q.split('|', 1) if '|' in q else ('', q)
        ss, se = int(d['sstart']), int(d['send'])
        hits.append({'accession': acc, 'taxon_group': grp,
                     'replicon': rest.split('_')[0], 'query_id': rest,
                     'pident': pid, 'qcovs': cov, 'length': int(d['length']),
                     'start': min(ss, se), 'end': max(ss, se),
                     'strand': '-' if ss > se else '+',
                     'bitscore': float(d['bitscore'])})

    # overlapping replicon hits on the same locus -> keep the higher identity
    hits.sort(key=lambda h: (-h['pident'], -h['bitscore']))
    kept = []
    for h in hits:
        if not any(min(h['end'], k['end']) - max(h['start'], k['start']) > 0
                   for k in kept):
            kept.append(h)
    kept.sort(key=lambda h: h['start'])

    cols = ['accession', 'replicon', 'taxon_group', 'query_id', 'pident',
            'qcovs', 'length', 'start', 'end', 'strand', 'bitscore']
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for h in kept:
            fh.write('\t'.join(str(h[c]) for c in cols) + '\n')
    names = ','.join(h['replicon'] for h in kept) or '-'
    print(f'{acc}: {len(kept)} replicon(s): {names}')


if __name__ == '__main__':
    main()
