#!/usr/bin/env python3
"""
Window extraction around a relBE locus -- the F2 step of PIPELINE.md.

Locator-blind: it does not care whether the locus was found by seed BLAST or by
GenBank locus tag. Callers pass coordinates; this module only cuts and reorients.

Replaces the windowing half of scripts/extract_tadb3_loci.py, fixing three bugs
that would otherwise be baked into the +/-20 kb regeneration:

  1. No circular wrap-around. The old clamp `max(1, lo-FLANK), min(len, hi+FLANK)`
     flagged 6 loci on large circular replicons as truncated when nothing was
     truncated at all -- NC_022996 is 124,908 bp and was flagged `truncated=left`
     purely because the locus sits near coordinate 2,523.
  2. Truncation recorded in the wrong coordinate frame. `trunc` was computed as
     left/right in REPLICON coordinates, while flank_up/flank_dn were swapped for
     minus-strand loci. The two disagreed for every minus-strand locus.
  3. No whole-replicon case. Replicons shorter than one window (NC_019238 at
     5,637 bp) must emit the entire molecule once, not be clamped twice.

PIPELINE.md: truncation is recorded, never used to exclude. These functions never
drop a locus; they describe it.
"""

_COMP = str.maketrans('ACGTRYKMSWBDHVNacgtrykmswbdhvn',
                      'TGCAYRMKSWVHDBNtgcayrmkswvhdbn')


def rc(s):
    """Reverse complement, IUPAC-aware."""
    return s.translate(_COMP)[::-1]


def extract_window(seq, lo, hi, relB_is_complement, flank,
                   topology='circular', replicon_len=None):
    """
    Cut a window of `flank` bp either side of the locus spanning [lo, hi]
    (1-based, inclusive, replicon coordinates) and reorient so relB is on +.

    Returns a dict. All flank/truncation fields are POST-reorientation, i.e. they
    describe the emitted sequence, not the replicon.

      seq_out            the window, relB on the + strand
      rel_start/rel_end  locus position within seq_out (1-based, inclusive)
      flank_up/flank_dn  bases actually obtained upstream/downstream
      truncated_up/dn    requested flank not achieved on that side
      wrapped            window crosses the replicon origin
      whole_replicon     window is the entire replicon
      strand             'minus' if reoriented, else 'plus'
      orig               human-readable source interval
    """
    n = replicon_len if replicon_len is not None else len(seq)
    if not 1 <= lo <= hi <= n:
        raise ValueError(f'locus {lo}-{hi} outside replicon of length {n}')
    span = hi - lo + 1
    circular = str(topology).lower().startswith('circ')

    wrapped = False
    whole_replicon = False

    if circular and span + 2 * flank >= n:
        # Replicon is smaller than the requested window. Emit it once, rotated so
        # the locus keeps as much symmetric context as the molecule allows.
        whole_replicon = True
        spare = n - span
        up = spare // 2
        start0 = (lo - 1 - up) % n
        win = seq[start0:] + seq[:start0]
        up_rep = up
        dn_rep = spare - up
        wrapped = start0 != 0
        orig = f'{1}-{n}(whole,rot@{start0 + 1})'

    elif circular:
        start0 = (lo - 1 - flank) % n
        end0 = (hi + flank) % n            # exclusive end, 0-based
        if start0 < end0:
            win = seq[start0:end0]
        else:
            win = seq[start0:] + seq[:end0]
            wrapped = True
        up_rep = dn_rep = flank            # circular: full flank always available
        orig = f'{start0 + 1}-{end0 if end0 else n}' + ('(wrap)' if wrapped else '')

    else:  # linear -- clamp, and truncation here is real
        ws = max(1, lo - flank)
        we = min(n, hi + flank)
        win = seq[ws - 1:we]
        up_rep = lo - ws
        dn_rep = we - hi
        orig = f'{ws}-{we}'

    # Truncation in replicon frame, then reorient everything together.
    trunc_up_rep = up_rep < flank
    trunc_dn_rep = dn_rep < flank

    if relB_is_complement:
        win = rc(win)
        up_f, dn_f = dn_rep, up_rep
        trunc_up, trunc_dn = trunc_dn_rep, trunc_up_rep
        strand = 'minus'
    else:
        up_f, dn_f = up_rep, dn_rep
        trunc_up, trunc_dn = trunc_up_rep, trunc_dn_rep
        strand = 'plus'

    return {
        'seq_out': win,
        'rel_start': up_f + 1,
        'rel_end': up_f + span,
        'flank_up': up_f,
        'flank_dn': dn_f,
        'truncated_up': trunc_up,
        'truncated_dn': trunc_dn,
        'wrapped': wrapped,
        'whole_replicon': whole_replicon,
        'strand': strand,
        'window_len': len(win),
        'orig': orig,
    }
