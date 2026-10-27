#!/usr/bin/env python3
"""
Resolve external executables so the same scripts run on a laptop and a cluster.

Resolution order, first hit wins:
  1. an explicit --blast-bin / bin_dir argument
  2. $BLAST_BIN
  3. PATH

Local runs set an explicit directory because BLAST+ here is a native arm64
build at ~/opt that is deliberately NOT on PATH (Homebrew on this machine is
x86_64/Rosetta, so brew's BLAST would be a downgrade).

Cluster runs leave it empty: the binary comes from an Apptainer container or a
module, and PATH is the right place to look. Hardcoding a macOS path makes the
pipeline unrunnable anywhere else.
"""

import os
import shutil


def resolve(name, bin_dir=''):
    """Return an absolute path to `name`, or raise with an actionable message."""
    for d in (bin_dir, os.environ.get('BLAST_BIN', '')):
        if d:
            p = os.path.join(d, name)
            if os.path.isfile(p) and os.access(p, os.X_OK):
                return p
    p = shutil.which(name)
    if p:
        return p
    raise SystemExit(
        f"{name!r} not found.\n"
        f"  looked in: --blast-bin={bin_dir or '(unset)'}, "
        f"$BLAST_BIN={os.environ.get('BLAST_BIN', '(unset)')}, and PATH\n"
        f"  local  : pass --blast-bin /path/to/ncbi-blast-X/bin\n"
        f"  cluster: load a module or use a container that provides {name}")
