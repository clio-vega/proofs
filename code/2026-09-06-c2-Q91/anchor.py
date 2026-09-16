"""Independent anchor: R_e(-1) = M_{p_e}, with p_e-multiplication computed by a
THIRD engine (symfunc.py: monomial basis + Kostka, no MN, no abacus)."""
import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-03-Q75')
import symfunc as S
from bead import R_bead, trim, t

ok = 0; tot = 0; bad = []
for e in (2, 3, 4):
    pe = S.p_m(e)
    for n in range(0, 7):
        for lam in (S.parts_of(n) if n else [()]):
            lam = trim(lam)
            nv = n + e + 1
            prod = S.mult(S.schur_m(lam), pe, nv)      # p_e * s_lam in monomials
            lhs = {k: v for k, v in S.to_schur(prod).items() if v != 0}
            rhs = {k: sp.expand(v.subs(t, -1)) for k, v in R_bead(lam, e).items()}
            rhs = {k: v for k, v in rhs.items() if v != 0}
            tot += 1
            if lhs == rhs: ok += 1
            else: bad.append((e, lam, lhs, rhs))
print(f"ANCHOR  R_e(-1) == M_(p_e)  [symfunc monomial engine]:  {ok}/{tot}")
for b in bad[:3]: print("  MISMATCH", b)
