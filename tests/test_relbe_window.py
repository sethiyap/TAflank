import sys, os, random
_R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), "..", "scripts"))
from relbe_window import extract_window, rc

random.seed(7)
def synth(n): return ''.join(random.choice('ACGT') for _ in range(n))

FAIL = []
def _abbr(v, n=24):
    r = repr(v)
    return r if len(r) <= n else r[:n-9] + '..' + r[-6:] + f'({len(v)}bp)'

def check(name, got, want):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<46} {_abbr(got):>26}")
    if not ok:
        print(f"          expected: {_abbr(want)}")
        FAIL.append(name)

print("=== 1. NC_022996-like: circular 124,908 bp, locus at 2523 (was FAKE truncated=left) ===")
n = 124908; lo = 2523; hi = lo + 526
r = extract_window(synth(n), lo, hi, False, 20000, 'circular', n)
check("wrapped", r['wrapped'], True)
check("truncated_up (was wrongly True)", r['truncated_up'], False)
check("truncated_dn", r['truncated_dn'], False)
check("flank_up", r['flank_up'], 20000)
check("window_len", r['window_len'], 527 + 40000)
check("whole_replicon", r['whole_replicon'], False)

print("\n=== 2. NC_019238-like: circular 5,637 bp -- smaller than one window ===")
n = 5637; lo = 2000; hi = lo + 526
r = extract_window(synth(n), lo, hi, False, 20000, 'circular', n)
check("whole_replicon", r['whole_replicon'], True)
check("window_len == replicon_len", r['window_len'], n)
check("truncated_up", r['truncated_up'], True)
check("no sequence duplicated", r['window_len'] <= n, True)
check("locus intact in window", r['rel_end'] - r['rel_start'] + 1, 527)

print("\n=== 3. KF033132-like: LINEAR 39,278 bp, minus strand, real truncation ===")
n = 39278; lo = 30000; hi = lo + 526   # near right end -> dn truncated in replicon frame
r = extract_window(synth(n), lo, hi, True, 20000, 'linear', n)
check("strand", r['strand'], 'minus')
check("replicon-right trunc -> truncated_UP after flip", r['truncated_up'], True)
check("truncated_dn", r['truncated_dn'], False)
check("flank_up (= replicon dn flank)", r['flank_up'], n - hi)
check("flank_dn (full)", r['flank_dn'], 20000)
check("wrapped (linear never wraps)", r['wrapped'], False)

print("\n=== 4. REAL CP018961: circular 289,903, relBE 257696-258222, minus, flank=20000 ===")
seq = ''.join(l.strip() for l in open(_R+'/data/public_16/fasta/CP018961.1.fna') if not l.startswith('>'))
print(f"  loaded {len(seq)} bp")
r = extract_window(seq, 257696, 258222, True, 20000, 'circular', len(seq))
check("window_len", r['window_len'], 40527)
check("relBE_in_window start", r['rel_start'], 20001)
check("relBE_in_window end", r['rel_end'], 20527)
check("flank_up", r['flank_up'], 20000)
check("flank_dn", r['flank_dn'], 20000)
check("not truncated", (r['truncated_up'], r['truncated_dn']), (False, False))
seed = ''.join(l.strip() for l in open(_R+'/anchors/relBE_seed.fa') if not l.startswith('>'))
check("seed sequence sits exactly at 20001-20527",
      r['seq_out'][r['rel_start']-1:r['rel_end']], seed)

print("\n=== 5. round-trip: plus-strand locus is unchanged by rc() ===")
n=50000; s=synth(n); lo=25000; hi=lo+526
a = extract_window(s, lo, hi, False, 20000, 'circular', n)
b = extract_window(s, lo, hi, True,  20000, 'circular', n)
check("minus window is rc of plus window", b['seq_out'], rc(a['seq_out']))

print("\n" + ("ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}"))
sys.exit(1 if FAIL else 0)
