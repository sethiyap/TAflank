#!/usr/bin/env python3
"""
Print named columns of a TSV, aligned.

  tsvcols.py results/cargo/CP028197*.cargo.tsv family source cov offset_from_bait
  tsvcols.py results/cargo/CP028197*.cargo.tsv            # list what is there

Exists because `cut -f2,8,9,10 | column -t` silently returns the WRONG columns
the moment a column is inserted ahead of one of those numbers -- which is
exactly what happened when clipped_at_window_edge landed at position 10 of the
cargo table and offset_from_bait slid out of the view. Three positional-index
bugs in the scripts themselves had the same cause. Asking for columns by name
cannot go quietly wrong: a name that is not there is an error, not a different
column.

Several files may be given; the header is printed once and rows are shown in
file order, so a glob over every locus gives one combined table.
"""

import sys


def read(path):
    with open(path) as fh:
        rows = [l.rstrip('\n').split('\t') for l in fh if l.strip()]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def main(argv):
    if not argv:
        sys.exit(__doc__.strip())

    files = [x for x in argv if '\t' not in x and
             (x.endswith('.tsv') or x.endswith('.txt') or '/' in x)]
    names = [x for x in argv if x not in files]
    if not files:
        sys.exit('no input file given')

    hdr, _ = read(files[0])
    if not hdr:
        sys.exit(f'{files[0]}: empty')

    if not names:
        print(f'{len(hdr)} columns in {files[0]}:')
        for i, c in enumerate(hdr, 1):
            print(f'  {i:>3}  {c}')
        return

    missing = [n for n in names if n not in hdr]
    if missing:
        sys.exit(f'no such column: {", ".join(missing)}\n'
                 f'available: {", ".join(hdr)}')

    out = [names]
    for f in files:
        h, rows = read(f)
        if h != hdr:
            bad = [n for n in names if n not in h]
            if bad:
                sys.exit(f'{f}: missing column(s) {", ".join(bad)}')
        ix = [h.index(n) for n in names]
        for r in rows:
            out.append([r[i] if i < len(r) else '' for i in ix])

    w = [max(len(r[i]) for r in out) for i in range(len(names))]
    for r in out:
        print('  '.join(v.ljust(w[i]) for i, v in enumerate(r)).rstrip())


if __name__ == '__main__':
    main(sys.argv[1:])
