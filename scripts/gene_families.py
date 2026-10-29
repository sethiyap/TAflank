#!/usr/bin/env python3
"""
Stage 3 -- gene families, so the figure's CDS track can be coloured.

  gene_families.py --gff run_bakta/bakta --loci run_bakta/loci \
                   --out results/gene_families.tsv

WHY THIS EXISTS
---------------
`gene_blocks.py` keys CDS colour on the product string. 127 of the 314
Bakta-named CDS are "hypothetical protein", so they collapse onto one colour
and large stretches of every locus are painted identically -- which reads as
homology that was never tested. Colour has to come from sequence, not from a
name, and that means clustering the proteins.

PROTEINS ARE TRANSLATED FROM THE GFF3, NOT READ FROM BAKTA'S .faa
-----------------------------------------------------------------
The BAKTA process in main.nf declares only `${prefix}.gff3` as output, so
publishDir copied only that; the .faa stayed in the Nextflow work directory,
which lived in $PBS_JOBFS and is gone. Re-running Bakta to recover a file we
can derive exactly is not worth two hours of queue: the GFF3 carries the CDS
coordinates, strand and phase that Pyrodigal called, and translating them
against the same window FASTA Bakta was given reproduces the proteins.

The translation is verified rather than assumed -- `--check` asserts that
relB and relE come back as the sequences recorded in CLAUDE.md.

(main.nf should still publish the .faa for future runs; this path stays as
the way to rebuild families from a run that did not.)

CLUSTERING ENGINE
-----------------
`mmseqs` where available -- it is what the sweep will need, and matches the
intent of Acman's Roary `-i 0.9 -s`. Where it is not, `blastp` all-vs-all plus
single-linkage gives the same answer at this size (314 proteins is nothing)
and needs no install. `--engine` forces one; the default tries mmseqs first.

FAMILY IDs ARE CONTENT-ADDRESSED
--------------------------------
A cluster index renumbers whenever the dataset grows, so colours would shuffle
between runs and no figure would be comparable with the last. The family id is
a hash of its representative sequence instead: a family keeps its identity, and
its colour, as loci are added.
"""

import argparse, hashlib, os, re, shutil, subprocess, sys
from collections import defaultdict

BASES = 'TCAG'
AAS = ('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
CODON = {a + b + c: AAS[i * 16 + j * 4 + k]
         for i, a in enumerate(BASES)
         for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

# from CLAUDE.md "The anchor" -- used by --check to prove the translation
KNOWN = {'relB': ('MGSINLRIDDEL', 'TLDEL', 79),
         'relE': ('MAYFLDFDERAL', 'VKRIL', 95)}


def rc(s):
    return s.translate(COMP)[::-1]


def translate(nt):
    aa = [CODON.get(nt[i:i + 3].upper(), 'X') for i in range(0, len(nt) - 2, 3)]
    if aa and aa[-1] == '*':
        aa.pop()
    # Prodigal/Pyrodigal CDS legitimately begin GTG or TTG; the initiator is
    # formyl-methionine regardless, and Bakta's own .faa writes M.
    if aa:
        aa[0] = 'M'
    return ''.join(aa)


def read_fasta_one(path):
    seq = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            if line.startswith('>'):
                break
            seq.append(line.strip())
    return ''.join(seq)


def parse_gff(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                if line.startswith('##FASTA'):
                    break
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            attr = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            m = re.search(r'UniRef:(UniRef50_\w+)', f[8])
            out.append({
                'gene_id': attr.get('ID', ''),
                'start': int(f[3]), 'end': int(f[4]), 'strand': f[6],
                'phase': int(f[7]) if f[7].isdigit() else 0,
                'product': attr.get('product', 'hypothetical protein'),
                'uniref50': m.group(1) if m else '',
            })
    return out


def cluster_mmseqs(faa, tmp, exe, min_id, cov):
    pref = os.path.join(tmp, 'clu')
    subprocess.run([exe, 'easy-cluster', faa, pref, os.path.join(tmp, 'mm'),
                    '--min-seq-id', str(min_id), '-c', str(cov),
                    '--cov-mode', '0', '-v', '1'],
                   check=True, capture_output=True, text=True)
    members = defaultdict(list)
    with open(pref + '_cluster.tsv') as fh:
        for line in fh:
            rep, mem = line.split()
            members[rep].append(mem)
    return members


def cluster_blastp(faa, tmp, bin_dir, min_id, cov, seqs):
    """All-vs-all blastp, then single linkage over edges passing both cuts."""
    db = os.path.join(tmp, 'db')
    subprocess.run([os.path.join(bin_dir, 'makeblastdb'), '-in', faa,
                    '-dbtype', 'prot', '-out', db],
                   check=True, capture_output=True, text=True)
    r = subprocess.run(
        [os.path.join(bin_dir, 'blastp'), '-query', faa, '-db', db,
         '-outfmt', '6 qseqid sseqid pident length qlen slen',
         '-evalue', '1e-5', '-max_target_seqs', '10000',
         '-num_threads', '4'],
        check=True, capture_output=True, text=True)

    parent = {k: k for k in seqs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for line in r.stdout.splitlines():
        q, s, pid, ln, ql, sl = line.split('\t')
        if q == s:
            continue
        pid, ln, ql, sl = float(pid), int(ln), int(ql), int(sl)
        # coverage of BOTH sequences, matching mmseqs --cov-mode 0, so a short
        # domain match cannot drag two unrelated proteins into one family
        if pid < min_id * 100 or ln < cov * ql or ln < cov * sl:
            continue
        a, b = find(q), find(s)
        if a != b:
            parent[a] = b

    members = defaultdict(list)
    for k in seqs:
        members[find(k)].append(k)
    return members


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gff', required=True, help='directory of Bakta *.gff3')
    p.add_argument('--loci', required=True, help='directory of window *.fa')
    p.add_argument('--out', required=True)
    p.add_argument('--faa', default='', help='also write the proteins here')
    p.add_argument('--engine', choices=['auto', 'mmseqs', 'blastp'],
                   default='auto')
    p.add_argument('--mmseqs', default=os.environ.get('MMSEQS', 'mmseqs'))
    p.add_argument('--blast-bin',
                   default=os.environ.get(
                       'BLAST_BIN',
                       os.path.expanduser('~/opt/ncbi-blast-2.17.0+/bin')))
    p.add_argument('--min-seq-id', type=float, default=0.90)
    p.add_argument('--cov', type=float, default=0.80)
    p.add_argument('--check', action='store_true',
                   help='assert relB/relE translate to the CLAUDE.md sequences')
    p.add_argument('--tmp', default='')
    a = p.parse_args()

    gffs = sorted(f for f in os.listdir(a.gff) if f.endswith('.gff3'))
    if not gffs:
        sys.exit(f'no *.gff3 under {a.gff}')

    genes, seqs = [], {}
    for g in gffs:
        lid = g[:-len('.gff3')]
        fa = os.path.join(a.loci, lid + '.fa')
        if not os.path.exists(fa):
            sys.exit(f'no window FASTA for {lid}: expected {fa}')
        win = read_fasta_one(fa)
        for c in parse_gff(os.path.join(a.gff, g)):
            nt = win[c['start'] - 1:c['end']]
            if c['strand'] == '-':
                nt = rc(nt)
            nt = nt[c['phase']:]
            aa = translate(nt)
            key = f"{lid}|{c['gene_id']}"
            c['locus_id'], c['key'], c['aa'] = lid, key, aa
            c['aa_len'] = len(aa)
            c['internal_stop'] = 'yes' if '*' in aa else 'no'
            genes.append(c)
            seqs[key] = aa

    nstop = sum(1 for c in genes if c['internal_stop'] == 'yes')
    if nstop:
        print(f'WARNING: {nstop} CDS translate with an internal stop codon; '
              'phase or window may be wrong', file=sys.stderr)

    if a.check:
        hits = {}
        for c in genes:
            for name, (head, tail, ln) in KNOWN.items():
                if c['aa'].startswith(head) and c['aa'].endswith(tail):
                    hits.setdefault(name, []).append(c['aa_len'])
        for name, (_, _, ln) in KNOWN.items():
            got = hits.get(name, [])
            assert got, f'--check: no CDS translated to {name}'
            assert ln in got, f'--check: {name} length {got} != {ln}'
            print(f'  check {name}: found {len(got)}, length {ln} aa  OK')

    tmp = a.tmp or os.path.join(os.path.dirname(a.out) or '.', '_families_tmp')
    os.makedirs(tmp, exist_ok=True)
    faa = a.faa or os.path.join(tmp, 'proteins.faa')
    os.makedirs(os.path.dirname(faa) or '.', exist_ok=True)
    with open(faa, 'w') as fh:
        for c in genes:
            fh.write(f">{c['key']}\n{c['aa']}\n")

    engine = a.engine
    exe = shutil.which(a.mmseqs) or (a.mmseqs if os.path.exists(a.mmseqs) else '')
    if engine == 'auto':
        engine = 'mmseqs' if exe else 'blastp'
    if engine == 'mmseqs':
        if not exe:
            sys.exit(f'--engine mmseqs but {a.mmseqs!r} not found')
        members = cluster_mmseqs(faa, tmp, exe, a.min_seq_id, a.cov)
    else:
        members = cluster_blastp(faa, tmp, a.blast_bin, a.min_seq_id, a.cov,
                                 seqs)

    # content-addressed id: hash of the representative sequence, so a family
    # keeps its identity (and colour) when the dataset grows
    fam_of, fam_size, fam_rep = {}, {}, {}
    for rep, mem in members.items():
        fid = 'F' + hashlib.sha1(seqs[rep].encode()).hexdigest()[:8]
        for m in mem:
            fam_of[m] = fid
        fam_size[fid] = len(mem)
        fam_rep[fid] = rep

    cols = ['locus_id', 'gene_id', 'start', 'end', 'strand', 'aa_len',
            'family_id', 'family_size', 'family_n_loci', 'is_representative',
            'product', 'uniref50', 'internal_stop']
    n_loci_in = defaultdict(set)
    for c in genes:
        n_loci_in[fam_of[c['key']]].add(c['locus_id'])

    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for c in genes:
            fid = fam_of[c['key']]
            fh.write('\t'.join(str(v) for v in [
                c['locus_id'], c['gene_id'], c['start'], c['end'], c['strand'],
                c['aa_len'], fid, fam_size[fid], len(n_loci_in[fid]),
                'yes' if fam_rep[fid] == c['key'] else 'no',
                c['product'], c['uniref50'], c['internal_stop']]) + '\n')

    shared = sum(1 for f in n_loci_in if len(n_loci_in[f]) > 1)
    hypo = sum(1 for c in genes if 'hypothetical' in c['product'].lower())
    hypo_fams = {fam_of[c['key']] for c in genes
                 if 'hypothetical' in c['product'].lower()}
    print(f'{len(gffs)} loci, {len(genes)} CDS -> {len(fam_size)} families '
          f'({engine}, id>={a.min_seq_id}, cov>={a.cov})')
    print(f'  {shared} families span more than one locus')
    print(f'  {hypo} "hypothetical protein" CDS -> {len(hypo_fams)} families '
          f'(was 1 colour)')
    print(f'-> {a.out}')


if __name__ == '__main__':
    main()
