"""REFUSAL PANEL for CLAIM 1 (Omega M-convex), after the first attempt produced a
   control that was no control.

WHAT WENT WRONG the first time: I "dropped the pair structure" by replacing
l_i <= x_i+e_i <= h_i with l_i <= x_i <= h_i.  That is a BOX intersected with a
HYPERPLANE in Z^{2m} -- which is M-convex.  So the control was a positive
instance of the very theorem, and it returned 1037/1037 green.  Recorded.

A REAL control must break LAMINARITY, which is the only hypothesis Claim 1 uses:
the constraint family {x_i,e_i} (pairs) + {singletons} + {everything} is laminar.
Add a CROSSING constraint and M-convexity must fail.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itertools import product
from collections import Counter
from verify_chain import is_Mconvex_pairs

def build(m, l, h, D, cross=None):
    """cross = (c_lo, c_hi): extra CROSSING bounds c_lo <= e_i + x_{i+1} <= c_hi."""
    out = []
    for x in product(*[range(0, max(h)+1) for _ in range(m)]):
        for e in product(*[range(0, max(h)+1) for _ in range(m)]):
            if sum(x)+sum(e) != D: continue
            if any(not (l[i] <= x[i]+e[i] <= h[i]) for i in range(m)): continue
            if cross is not None:
                clo, chi = cross
                if any(not (clo <= e[i] + x[(i+1) % m] <= chi) for i in range(m)): continue
            out.append((x, e))
    return out

for name, cross in (("laminar (Claim 1 itself) -- must NOT refuse", None),
                    ("CROSSING bounds 0<=e_i+x_{i+1}<=1 -- must refuse", (0, 1)),
                    ("CROSSING bounds 1<=e_i+x_{i+1}<=2 -- must refuse", (1, 2)),
                    ("CROSSING bounds 0<=e_i+x_{i+1}<=2 -- must refuse", (0, 2))):
    st = Counter()
    for m in (2, 3):
        for l in product(range(0, 3), repeat=m):
            for h in product(*[range(l[i], 4) for i in range(m)]):
                for D in range(sum(l), sum(h)+1):
                    S = build(m, list(l), list(h), D, cross)
                    if len(S) < 2 or len(S) > 120: continue
                    st['Mconvex' if is_Mconvex_pairs(S) is None else 'REFUSED'] += 1
    print(f"  {name:52s} {dict(st)}", flush=True)
