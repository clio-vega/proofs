"""
Job A, c=4,5 confirmation via the exact 3-row Jacobi-Trudi determinant (dj.Mj_jt).
No clean pure-Python interior closed form for c=4,5, so we use the determinant at
modest m to confirm: above the box-fits gate (j0+2c<=b) the boundary indices
j in {b+1,...,b+c} strictly lose the 2-adic minimum, and the only boundary ties
are a finite small-b list.
"""
import sys
sys.path.insert(0, '/home/clio/projects/code/threerow-c1')
from dj import Mj_jt
from mn import Mj
from math import comb

def v2(n):
    n = abs(int(n))
    if n == 0:
        return float('inf')
    c = 0
    while n % 2 == 0:
        n //= 2; c += 1
    return c

def val(lam, j, m):
    M = Mj_jt(lam, j, m)
    if M == 0 or j > m:
        return None
    return j + 2 * (v2(comb(m, j)) + v2(M))

# Mj_jt (Jacobi-Trudi determinant) is verified == MN on import of dj.py
# ("D_j/JT formula check, mismatches = 0"), so we use it directly here.

for c in (4, 5):
    in_ties = 0
    in_min_delta = None
    in_min_at = None
    genuine_tie_b = set()
    shapes = 0
    mmax = 40 if c == 4 else 34
    for m in range(c, mmax + 1):
        n = 2 * m
        for a in range(1, n + 1):
            b = n - a - c
            if b < c or b > a:
                continue
            shapes += 1
            vals = {}
            for j in range(0, b + c + 1):
                v = val((a, b, c), j, m)
                if v is not None:
                    vals[j] = v
            if not vals:
                continue
            V = min(vals.values())
            Jstar = [j for j, v in vals.items() if v == V]
            j0 = min(jj for jj in vals if jj <= b and vals[jj] == V) if any(jj <= b and vals[jj] == V for jj in vals) else min(Jstar)
            # boundary min
            bvals = [vals[j] for j in vals if j > b]
            interior_min = min((vals[j] for j in vals if j <= b), default=None)
            if any(j > b for j in Jstar):
                genuine_tie_b.add(b)
            if interior_min is not None and bvals:
                if j0 + 2 * c <= b:  # in-gate
                    d = min(bvals) - interior_min
                    if in_min_delta is None or d < in_min_delta:
                        in_min_delta = d; in_min_at = (a, b, c, m)
                    if d <= 0:
                        in_ties += 1
    print(f"c={c}: {shapes} shapes m<={mmax};  IN-GATE min Delta={in_min_delta} at {in_min_at}, "
          f"in-gate ties={in_ties};  genuine boundary-tie b-values={sorted(genuine_tie_b)}")
