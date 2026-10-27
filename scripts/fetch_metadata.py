#!/usr/bin/env python3
"""
Tier B -- fetch record metadata and classify replicon type.

  fetch_metadata.py --accessions results/bait_hits.tsv --col sacc \
                    --out results/record_metadata.tsv [--offline]

Split into a networked fetch and an offline build, deliberately: cluster compute
nodes may have no outbound access (Gadi's standard queues do not), so the fetch
populates a versioned on-disk cache and --offline rebuilds the table from it
with no network at all. Cache the cache, not the answer.

Classification uses NCBI's STRUCTURED esummary fields, not title parsing:
  genome        -> plasmid | chromosome | ''      (the replicon call)
  completeness  -> complete | partial
  slen          -> record length

Records where `genome` is empty get replicon_type=unknown, NOT forced into a
bucket. In a broad sweep many hits land on contigs from draft assemblies where
NCBI genuinely cannot say. PlasmidFinder provides independent evidence for
those -- a replicon hit says plasmid even when the metadata is silent.
"""

import argparse, json, os, sys, time, urllib.parse, urllib.request

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
FIELDS = ['accessionversion', 'title', 'organism', 'taxid', 'slen', 'genome',
          'completeness', 'biomol', 'moltype', 'topology', 'strain',
          'sourcedb', 'createdate', 'updatedate', 'assemblygi', 'biosample']


def api_key():
    p = os.path.expanduser('~/.ncbi_api_key')
    return open(p).read().strip() if os.path.exists(p) else ''


def esummary(accs, delay):
    key = api_key()
    q = {'db': 'nuccore', 'id': ','.join(accs), 'retmode': 'json'}
    if key:
        q['api_key'] = key
    url = EUTILS + '/esummary.fcgi?' + urllib.parse.urlencode(q)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)['result']
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
        finally:
            time.sleep(delay)
    return {}


def classify_replicon(r):
    """Return (replicon_type, replicon_name, evidence).

    The structured esummary `genome` field is authoritative where present, but
    it is NOT always populated: CP021960 is a 5.6 Mb circular K. pneumoniae
    "complete genome" with genome=None. Fall back to the submitter's title,
    which is explicit text, and record which route was taken so a reader can
    tell a structured call from an inferred one.

    Deliberately no size heuristic -- megaplasmids exceed 1 Mb and would be
    misassigned as chromosomes.
    """
    title = r.get('title', '') or ''
    low = title.lower()
    g = (r.get('genome') or '').lower()

    if g in ('plasmid', 'chromosome'):
        rtype, evid = g, 'esummary_genome'
    elif 'plasmid' in low:
        rtype, evid = 'plasmid', 'title'
    elif 'chromosome' in low or 'complete genome' in low:
        rtype, evid = 'chromosome', 'title'
    else:
        return 'unknown', '', 'none'

    if rtype == 'plasmid':
        rname = title.split('plasmid', 1)[1] if 'plasmid' in low else ''
        rname = rname.strip(' :,.').split(',')[0].strip()
    else:
        rname = 'chromosome'
    return rtype, rname, evid


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--accessions', required=True)
    p.add_argument('--col', default='sacc')
    p.add_argument('--out', required=True)
    p.add_argument('--cache', default='data/metadata_cache')
    p.add_argument('--batch', type=int, default=200)
    p.add_argument('--offline', action='store_true',
                   help='build from cache only; never touch the network')
    a = p.parse_args()

    with open(a.accessions) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        if a.col not in hdr:
            sys.exit(f'column {a.col!r} not in {a.accessions}; have {hdr}')
        i = hdr.index(a.col)
        accs = sorted({l.split('\t')[i] for l in fh if l.strip()})

    os.makedirs(a.cache, exist_ok=True)
    need = [x for x in accs
            if not os.path.exists(os.path.join(a.cache, x + '.json'))]

    if need and a.offline:
        sys.exit(f'--offline but {len(need)} accessions are not cached: '
                 f'{need[:5]}... run once without --offline')

    delay = 0.12 if api_key() else 0.35     # 10/s with a key, 3/s without
    for b in range(0, len(need), a.batch):
        chunk = need[b:b + a.batch]
        print(f'  fetching {b + 1}-{b + len(chunk)} of {len(need)}', flush=True)
        res = esummary(chunk, delay)
        for uid in res.get('uids', []):
            rec = res[uid]
            acc = rec.get('accessionversion', '')
            for name in {acc, acc.split('.')[0]} & set(chunk) or {acc}:
                with open(os.path.join(a.cache, name + '.json'), 'w') as fh:
                    json.dump(rec, fh)

    out_hdr = ['accession'] + FIELDS + ['replicon_type', 'replicon_name',
                                        'genus', 'is_complete',
                                        'replicon_type_evidence']
    rows = []
    for x in accs:
        f = os.path.join(a.cache, x + '.json')
        if not os.path.exists(f):
            rows.append([x] + [''] * len(FIELDS) + ['unknown', '', '', '', 'no_record'])
            continue
        r = json.load(open(f))
        vals = [r.get(k, '') for k in FIELDS]
        rtype, rname, evid = classify_replicon(r)
        org = r.get('organism', '')
        rows.append([x] + vals + [rtype, rname, org.split(' ')[0] if org else '',
                                  'yes' if (r.get('completeness') or '').lower()
                                  == 'complete' else 'no', evid])

    with open(a.out, 'w') as fh:
        fh.write('\t'.join(out_hdr) + '\n')
        for r in rows:
            fh.write('\t'.join(str(v) for v in r) + '\n')

    n = {}
    ti = out_hdr.index('replicon_type')
    for r in rows:
        n[r[ti]] = n.get(r[ti], 0) + 1
    print(f'{len(rows)} records -> {a.out}')
    for k in sorted(n):
        print(f'  {k:<12} {n[k]}')


if __name__ == '__main__':
    main()
