#!/usr/bin/env python3
"""
Run pangraph on the SAME paralog/completeness-filtered loci windows the
blastn pipeline already built -- this does not re-do that filtering, it
consumes its output (run_all_v2/loci, or wherever --loci points).

  build_pangraph_graph.py --loci run_all_v2/loci --out results/pangraph/graph.json

pangraph is a true N-way alignment graph (Noll, Molari & Neher 2023,
Microbial Genomics 9:001034; Rust rewrite, github.com/neherlab/pangraph) --
every input genome is an explicit ordered walk over blocks, each block with
its own strand, so a rearrangement is represented directly rather than
approximated. The blastn-based paralog/completeness filter in
bait_search.py stays exactly as is -- this script only changes how
homology BLOCKS are found among the loci that filter already kept.

Separate from scripts/pangraph_blocks.py on purpose, same split as
bait_search.py vs select_loci.py elsewhere in this repo: one script runs
the external tool and writes a flat file, a different one reads that file
and draws conclusions from it.
"""
import argparse, glob, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolpath import resolve  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--loci', required=True, help='dir of locus window *.fa')
    p.add_argument('--out', required=True, help='output graph.json path')
    p.add_argument('--workdir', default='tmp/pangraph',
                   help='combined input fasta lives here -- NOT /tmp, per '
                        'configs/system_paths.yaml temp_dirs.temp_root')
    p.add_argument('--jobs', type=int, default=0,
                   help='0 = pangraph default (all available CPU threads)')
    p.add_argument('--pangraph-bin', default=os.environ.get('PANGRAPH_BIN', ''),
                   help='dir holding the pangraph binary; empty => PATH')
    a = p.parse_args()

    files = sorted(glob.glob(os.path.join(a.loci, '*.fa')))
    if not files:
        sys.exit(f'no *.fa under {a.loci}')
    os.makedirs(a.workdir, exist_ok=True)
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)

    combined = os.path.join(a.workdir, 'combined.fa')
    with open(combined, 'w') as out:
        for f in files:
            with open(f) as fh:
                out.write(fh.read())

    cmd = [resolve('pangraph', a.pangraph_bin), 'build']
    if a.jobs:
        cmd += ['-j', str(a.jobs)]
    cmd += [combined, '-o', a.out]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        sys.exit('pangraph build failed:\n' + res.stderr)

    print(f'{len(files)} loci -> {a.out}')
    if res.stderr.strip():
        print(res.stderr.strip())


if __name__ == '__main__':
    main()
