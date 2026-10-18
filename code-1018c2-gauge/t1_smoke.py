"""Smoke test: reproduce (SRC) d=1 and (QC) d=2 through the generic AtilU route."""
import sympy as sp, time
from gauge_check import *

t0 = time.time()
x = sp.symbols('x')
N = 2
w = wvars(N)

# ---- (SRC), type C: alphabet w_2 = eps, D = [eps], R = [w_1], I = () and (1,)
print("=== d=1 boundary (SRC), N=2, type C ===")
for eps in (1, -1):
    wl = [w[0], sp.Integer(eps)]
    for I in [(), (1,), (2,)]:
        res = adjunction(N, [w[0]], [sp.Integer(eps)], wl, I, 'C', 2)
        print(f"  eps={eps:+d} I={I}: {res}")

# ---- (QC), type C: alphabet w = (x, 1/x), D=[x,1/x], R=[], I=()
print("=== d=2 pair (QC), N=2, type C, R empty ===")
wl = [x, 1/x]
for I in [(), (1,), (2,)]:
    res = adjunction(N, [], [x, 1/x], wl, I, 'C', 2)
    print(f"  I={I}: {res}")
print(f"[{time.time()-t0:.1f}s]")
