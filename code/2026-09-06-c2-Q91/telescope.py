"""Verify the telescoping potential in the proof of [R_e^*,R_e]=[e]_{t^2}.

  f(b) := (m(b)-m(b+e)) t^{2 N_b},     N_b = #(M cap (b,b+e))
  g(b) := sum_{i=0}^{e-1} m(b+i) t^{2 #(M cap (b+i, b+e))}
  CLAIM: f(b) = g(b) - g(b+1) for every b, and g(-inf)=[e]_{t^2}, g(+inf)=0.
"""
import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine
from bead import maya, trim, t

def check(lam, e):
    lo = -(len(trim(lam)) + 3 * e + 8); M = maya(lam, lo)
    hi = (lam[0] if lam else 0) + 3 * e + 8
    m = lambda x: 1 if x in M else 0
    N = lambda b: sum(m(j) for j in range(b + 1, b + e))
    f = lambda b: (m(b) - m(b + e)) * t ** (2 * N(b))
    g = lambda b: sum(m(b + i) * t ** (2 * sum(m(j) for j in range(b + i + 1, b + e)))
                      for i in range(e))
    bad = 0
    for b in range(lo + e + 2, hi - e - 2):
        if sp.expand(f(b) - (g(b) - g(b + 1))) != 0:
            bad += 1
    S = sp.expand(sum(f(b) for b in range(lo + e + 2, hi - e - 2)))
    return bad, S, sp.expand(g(lo + e + 2)), sp.expand(g(hi - e - 2))

pred = lambda e: sp.expand(sum(t ** (2 * h) for h in range(e)))
tot = 0; bad_tot = 0; sbad = 0
for e in (1, 2, 3, 4, 5):
    for n in range(0, 8):
        for lam in (engine.parts_of(n) if n else ((),)):
            b, S, gm, gp = check(lam, e)
            tot += 1; bad_tot += b
            if S != pred(e) or gm != pred(e) or gp != 0: sbad += 1
print(f"telescoping identity f(b)=g(b)-g(b+1): {tot} (lam,e) pairs, "
      f"{bad_tot} failing b's")
print(f"sum f = g(-inf) = [e]_(t^2) and g(+inf)=0 : {tot-sbad}/{tot}")
