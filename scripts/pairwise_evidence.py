#!/usr/bin/env python3
"""
Stage 6 input -- reconcile the two independent lines of evidence into one
call per locus pair, so "ancestral vehicle vs. independent acquisition"
(CLAUDE.md Stage 6) is read off a table instead of eyeballed across two
figures.

  pairwise_evidence.py --genes results/structural_pairs.tsv \
                        --seq results/pangraph/pangraph_pairwise.tsv \
                        --out results/pairwise_evidence.tsv

structural_pairs.tsv (structural_groups.py) says whether two loci carry the
SAME NAMED elements (core_jaccard over the IS/cargo gene set: ISEsa1, czc,
relB, ...). pangraph_pairwise.tsv (pangraph_blocks.py) says whether two
loci's RAW SEQUENCE is actually alignable, independent of any name. Agreeing is
itself the finding:

  genes agree  + sequence agrees    -> CONCORDANT     same cargo, same DNA:
                                        real shared ancestry, one vehicle
  genes agree  + sequence disagrees -> DISCORDANT      same cargo, different
                                        DNA: convergent / independent
                                        acquisition of similar elements
  genes disagree + sequence agrees  -> SEQUENCE_ONLY   shared backbone the
                                        gene-based pass missed entirely
                                        (unnamed CDS, intergenic) -- a
                                        candidate for extending the anchor
  genes disagree + sequence disagrees -> UNRELATED     no evidence of
                                        shared history either way

Thresholds are arguments, not constants, because "agree" is a judgement call
that should be visible and re-runnable, not buried in the code.
"""
import argparse, csv, sys


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--genes', required=True,
                   help='results/structural_pairs.tsv (structural_groups.py)')
    p.add_argument('--seq', required=True,
                   help='results/pangraph/pangraph_pairwise.tsv (pangraph_blocks.py)')
    p.add_argument('--out', default='results/pairwise_evidence.tsv')
    p.add_argument('--gene-jaccard-min', type=float, default=0.8,
                   help='>= this core_jaccard counts as "genes agree"')
    p.add_argument('--seq-frac-min', type=float, default=0.10,
                   help='>= this fraction of the 2x flank shared as '
                        'contiguous sequence from the core counts as '
                        '"sequence agrees" -- deliberately low: even a '
                        'short but definite shared run right at the core '
                        'edge is evidence, this is not asking for the '
                        'whole window to match')
    a = p.parse_args()

    genes = {(r['locus_a'], r['locus_b']): float(r['core_jaccard'])
             for r in read_tsv(a.genes)}
    seq = {(r['locus_a'], r['locus_b']): r for r in read_tsv(a.seq)}

    pairs = sorted(set(genes) | set(seq))
    if set(genes) != set(seq):
        only_g = set(genes) - set(seq)
        only_s = set(seq) - set(genes)
        print(f'WARNING: {len(only_g)} pair(s) only in --genes, '
              f'{len(only_s)} only in --seq -- these two tables were not '
              f'built from the same locus set', file=sys.stderr)

    rows = []
    for key in pairs:
        gj = genes.get(key)
        sr = seq.get(key)
        sfrac = float(sr['shared_frac']) if sr else 0.0
        genes_agree = gj is not None and gj >= a.gene_jaccard_min
        seq_agrees = sfrac >= a.seq_frac_min
        if genes_agree and seq_agrees:
            call = 'CONCORDANT'
        elif genes_agree and not seq_agrees:
            call = 'DISCORDANT'
        elif not genes_agree and seq_agrees:
            call = 'SEQUENCE_ONLY'
        else:
            call = 'UNRELATED'
        rows.append({
            'locus_a': key[0], 'locus_b': key[1],
            'organism_short_a': sr.get('organism_short_a', '') if sr else '',
            'organism_short_b': sr.get('organism_short_b', '') if sr else '',
            'gene_jaccard': f'{gj:.3f}' if gj is not None else 'NA',
            'shared_up_bp': sr['shared_up_bp'] if sr else '0',
            'shared_down_bp': sr['shared_down_bp'] if sr else '0',
            'shared_frac': f'{sfrac:.4f}',
            'call': call,
        })

    with open(a.out, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()),
                           delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    counts = {}
    for r in rows:
        counts[r['call']] = counts.get(r['call'], 0) + 1
    print(f'{len(rows)} pairs -> {a.out}')
    for call in ('CONCORDANT', 'DISCORDANT', 'SEQUENCE_ONLY', 'UNRELATED'):
        n = counts.get(call, 0)
        print(f'  {call:<14} {n}')
        if call in ('CONCORDANT', 'DISCORDANT') and n:
            for r in rows:
                if r['call'] == call:
                    a_acc = r['locus_a'].split('__')[0]
                    b_acc = r['locus_b'].split('__')[0]
                    print(f'      {a_acc} - {b_acc}  '
                          f'genes={r["gene_jaccard"]}  '
                          f'seq_shared={r["shared_frac"]}')

    if counts.get('DISCORDANT'):
        print(f'\n{counts["DISCORDANT"]} pair(s) carry the same named '
              'elements on DIFFERENT DNA -- candidate evidence for '
              'independent/convergent acquisition rather than one '
              'ancestral vehicle (CLAUDE.md Stage 6 decisive test)')


if __name__ == '__main__':
    main()
