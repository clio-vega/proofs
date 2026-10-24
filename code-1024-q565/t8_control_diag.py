"""T8.  Diagnose the 19 zeros the 'hygienic' sampler produced.

The brief said a zero under Q559 hygiene = engine regression.  But the ORIGINAL
sampler (code-1023-q559/sampler.py) draws F(randint(2,40), randint(1,7)): every
w_a is POSITIVE.  Positivity is a FOURTH exclusion, unnamed in the hygiene note
and unnamed in today's brief.  Classify each zero and then re-run the control in
the positive orthant, where it must stay green.
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, B
from sampler import clean_w

def mech(w, r, theta):
    """classify a zero: M1 structural / M2 alternant / M3 cancellation."""
    n = len(w); m = [theta[j]+j for j in range(r)]
    Bs = {T: B(w, T) for T in itertools.combinations(range(n), r)}
    if all(v == 0 for v in Bs.values()): return "M1 structural (all B=0)"
    from psi import alternant_ratio
    alts = {}
    for T in Bs:
        x = [w[i] for i in T]
        alts[T] = alternant_ratio(x, m)
    live = [T for T in Bs if Bs[T] != 0]
    if all(alts[T] == 0 for T in live): return "M2 alternant (every live term has det=0)"
    return "M3 cancellation (nonzero terms summing to 0)"

random.seed(99)
VALS = [F(2),F(3),F(5),F(7),F(-2),F(-3),F(4),F(-5),F(11),F(-7),F(3,2),F(-5,3)]
zeros = []
for r in range(1,4):
    for n in range(r, r+3):
        for _ in range(60):
            w = random.sample(VALS, n)
            if any(x*x == 1 for x in w): continue
            if any(w[i]*w[j] == 1 for i in range(n) for j in range(i+1,n)): continue
            if len(set(w)) != n: continue
            theta = list(range(r))
            try: v = psi_headsum(w, theta)
            except ZeroDivisionError: continue
            if v == 0:
                assert psi_alt(w, theta) == 0, "implementations disagree!"
                zeros.append((n, r, tuple(theta), tuple(w)))

print(f"zeros to diagnose: {len(zeros)}")
from collections import Counter
cnt = Counter(); shown = set()
for (n, r, th, w) in zeros:
    k = mech(list(w), r, list(th))
    cnt[(n, r, k)] += 1
    key = (n, r, k)
    if key not in shown:
        shown.add(key)
        neg = [x for x in w if x < 0]
        e = F(1)
        for x in w: e *= x
        print(f"  n={n} r={r} theta={th}  {k}")
        print(f"     M={list(w)}   negatives={neg}   e={e}")
print()
print("classification totals:")
for (n, r, k), c in sorted(cnt.items()):
    print(f"  n={n} r={r:>1}  {k:42} x{c}")
print()
allneg = all(any(x < 0 for x in w) for (_,_,_,w) in zeros)
print(f"EVERY diagnosed zero contains a NEGATIVE w_a: {allneg}")

print()
print("=== the control, re-run in the POSITIVE orthant (the original sampler) ===")
rng = random.Random(2024)
zer = tot = 0
for r in range(1,4):
    for n in range(r, r+3):
        for _ in range(200):
            w = clean_w(n, rng)
            theta = list(range(r))
            try: v = psi_headsum(w, theta)
            except ZeroDivisionError: continue
            tot += 1; zer += (v == 0)
print(f"    positive hygienic sampler: {zer} zeros in {tot} draws  (must be 0)")
print()
print("=== and the SAME sampler with signs allowed: the suppressed mechanism ===")
rng2 = random.Random(2024)
zer2 = tot2 = 0; byr = Counter()
for r in range(1,4):
    for n in range(r, r+3):
        for _ in range(200):
            w = clean_w(n, rng2)
            w = [x if rng2.random() < 0.5 else -x for x in w]
            if len(set(w)) != n: continue
            if any(x*x == 1 for x in w): continue
            if any(w[i]*w[j] == 1 for i in range(n) for j in range(i+1,n)): continue
            theta = list(range(r))
            try: v = psi_headsum(w, theta)
            except ZeroDivisionError: continue
            tot2 += 1
            if v == 0:
                zer2 += 1; byr[(n, r, mech(w, r, theta))] += 1
print(f"    sign-randomised hygienic sampler: {zer2} zeros in {tot2} draws")
for k, c in sorted(byr.items()): print(f"      n={k[0]} r={k[1]}  {k[2]} x{c}")
