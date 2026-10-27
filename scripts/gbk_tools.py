#!/usr/bin/env python3
"""
Minimal GenBank reader -- stdlib only.

biopython is installed under no interpreter on this machine, so every script
here stays stdlib-only (as scripts/extract_tadb3_loci.py already did). This
module reads only what the pipeline needs from the flatfile header and the
`source` feature; it deliberately does not parse CDS features, because loci in
this stage are located by seed homology, not by locus tag.
"""

import re

# /country was renamed /geo_loc_name in 2023. Records in data/public_16 use the
# new spelling exclusively (12/16 have geo_loc_name, 0/16 have country), but the
# TADB3 parents are older, so both must be accepted.
_SOURCE_QUALS = ('organism', 'strain', 'plasmid', 'host', 'isolation_source',
                 'geo_loc_name', 'country', 'collection_date', 'lat_lon',
                 'serovar', 'isolate', 'db_xref', 'mol_type', 'note')


def parse_locus(path):
    """Read the LOCUS line and the source feature. Returns a dict."""
    acc, length, topology = None, None, None
    quals, in_source, buf, key = {}, False, [], None

    def flush():
        if key and buf:
            quals.setdefault(key, []).append(' '.join(buf).strip('"'))

    with open(path) as fh:
        for line in fh:
            if acc is None and line.startswith('LOCUS'):
                p = line.split()
                acc, length = p[1], int(p[2])
                topology = 'circular' if 'circular' in line else 'linear'
                continue
            if line.startswith('ORIGIN') or line.startswith('//'):
                break
            if re.match(r'^ {5}\S', line):          # a feature key at col 6
                flush(); buf, key = [], None
                in_source = line[5:].split()[0] == 'source'
                continue
            if in_source and re.match(r'^ {21}/', line):
                flush(); buf = []
                m = re.match(r'^ {21}/([^=]+)=?(.*)', line.rstrip('\n'))
                key = m.group(1)
                buf = [m.group(2)] if m.group(2) else []
                continue
            if in_source and key and re.match(r'^ {21}\S', line):
                buf.append(line.strip())
    flush()

    out = {'accession': acc, 'replicon_len': length, 'topology': topology}
    for q in _SOURCE_QUALS:
        v = quals.get(q, [])
        out[q] = v[0].strip('"') if v else ''
    # /plasmid present => plasmid; its absence in a chromosomal record is how
    # results/nine_loci_classification_METHODS.md assigned replicon_type.
    out['replicon_type'] = 'plasmid' if out['plasmid'] else 'chromosome'
    out['replicon_name'] = out['plasmid'] or 'chromosome'
    for d in quals.get('db_xref', []):
        if d.startswith('taxon:'):
            out['taxid'] = d.split(':', 1)[1]
    out.setdefault('taxid', '')
    g = out['organism'].split()
    out['genus'] = g[0] if g else ''
    out['species'] = ' '.join(g[:2]) if len(g) > 1 else out['organism']
    out['country'] = (out['geo_loc_name'] or out['country']).split(':')[0].strip()
    return out


def read_fasta(path):
    """Return (header, sequence) of the first record. Sequence uppercased."""
    hdr, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:
                    break
                hdr = line[1:].strip()
            else:
                seq.append(line.strip())
    return hdr, ''.join(seq).upper()


def write_tsv(path, header, rows):
    """Every emitted table gets a header -- results/stage1_is_hits.tsv is
    currently headerless raw BLAST and any consumer has to guess."""
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in rows:
            fh.write('\t'.join('' if v is None else str(v) for v in r) + '\n')
