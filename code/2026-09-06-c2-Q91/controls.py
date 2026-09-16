"""Planted-defect controls.  EACH MUST FAIL.  A control that passes means the
harness is measuring something other than what I think."""
import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
from phase0 import compare

CTRL = [
    ("C0  conjecture as stated      (open, N,   (-t)^N)", dict()),
    ("C1  CLOSED interval [b,b+e]   (closed, N, (-t)^N)", dict(interval='closed')),
    ("C2  exponent N+1              (open, N+1, (-t)^N)", dict(shift=1)),
    ("C3  sign flip  t^N            (open, N,   (+t)^N)", dict(tsign=1)),
]
for name, kw in CTRL:
    a, tot, bad = compare([2, 3, 4], 7, **kw)
    verdict = "AGREES" if a == tot else "DISAGREES"
    print(f"{name}:  {a}/{tot}  -> {verdict}")
    if bad:
        e, lam, lhs, rhs = bad[0]
        print(f"      smallest witness: e={e}, lam={lam}")
        for mu in sorted(set(lhs) | set(rhs)):
            if sp.expand(lhs.get(mu, 0) - rhs.get(mu, 0)) != 0:
                print(f"        mu={mu}: shape {lhs.get(mu,0)}   ctrl {rhs.get(mu,0)}")
