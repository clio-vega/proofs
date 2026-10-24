"""T11.  Exact threshold grid for the negation-pair vanishing at theta=(0,..,r-1).
For each (lambda, f) with n = 2*lambda + f, record for which r  Psi^C = 0 on all draws.
Also separate M2 (all terms dead, predicted r > n-lambda) from M3 (cancellation)."""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, B, alternant_ratio

random.seed(1111)
POOL = [F(2),F(3),F(5),F(7),F(11),F(13),F(4),F(9),F(8),F(5,2),F(7,3),F(9,4),F(6),F(10),F(17),F(19)]

def draw(lam, f):
    vals = random.sample(POOL, lam + f)
    w = []
    for v in vals[:lam]: w += [v, -v]
    w += vals[lam:]
    n = len(w)
    if len(set(w)) != n or any(x == 0 for x in w): return None
    if any(w[i]*w[j] == 1 for i in range(n) for j in range(i, n)): return None
    return w

print("rows: (lambda,f,n).  cell = r: Z-all/Z-M2/Z-M3 over 12 draws; '.' = no zeros")
print(f"{'lam':>3} {'f':>2} {'n':>2} | " + " ".join(f"r={r}" for r in range(1, 9)))
for lam in range(0, 4):
    for f in range(0, 4):
        n = 2*lam + f
        if n < 1 or n > 8: continue
        cells = []
        for r in range(1, 9):
            if r > n: cells.append("   -"); continue
            theta = list(range(r)); m = [theta[j]+j for j in range(r)]
            if len(set(m)) != r: cells.append("   ?"); continue
            z = t = m2 = m3 = 0
            for _ in range(12):
                w = draw(lam, f)
                if w is None: continue
                try: v = psi_headsum(w, theta)
                except ZeroDivisionError: continue
                t += 1
                if v != 0: continue
                z += 1
                live = []
                for T in itertools.combinations(range(n), r):
                    x = [w[i] for i in T]
                    if B(w, T)*alternant_ratio(x, m) != 0: live.append(T)
                if live: m3 += 1
                else: m2 += 1
            if t == 0: cells.append("   .")
            elif z == 0: cells.append("   .")
            elif z == t: cells.append(f"{'M2' if m2==z else ('M3' if m3==z else 'MX')}*")
            else: cells.append(f"{z}/{t}")
        print(f"{lam:>3} {f:>2} {n:>2} | " + " ".join(f"{c:>4}" for c in cells))
print()
print("legend: M2* = zero on every draw, all terms dead (alternant covering law)")
print("        M3* = zero on every draw, live terms cancel (NEW, not explained by any")
print("              covering law: a genuine cancellation IDENTITY)")
print("        .   = no zeros")
print()
print("alternant covering law predicts M2 exactly when r > n - lambda:")
for lam in range(0, 4):
    for f in range(0, 4):
        n = 2*lam + f
        if 1 <= n <= 8:
            thr = n - lam
            print(f"   lam={lam} f={f} n={n}: all-terms-dead predicted for r >= {thr+1}")
