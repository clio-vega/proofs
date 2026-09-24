import sys; sys.path.insert(0,'.')
from anTL import *
import sympy as sp

# (G3): at k=1 and k=n-1, is R_e(t) t-dependent at all?
print("k=1 and k=n-1: degree of R_e(t) in t, entrywise")
for n in range(3,9):
    for k in (1, n-1):
        for e in range(1,n):
            M = cyclic_ribbon_adder(e,n,k)
            degs=set()
            for S,row in M.items():
                for T,c in row.items():
                    p = sp.Poly(sp.expand(c), t)
                    degs.add((p.degree(), len(p.all_coeffs())-1-p.all_coeffs()[::-1].index(next(x for x in p.all_coeffs()[::-1] if x!=0))))
                    # simpler: record monomials present
            mons=set()
            for S,row in M.items():
                for T,c in row.items():
                    for mono in sp.Add.make_args(sp.expand(c)):
                        mons.add(sp.degree(mono, t))
            print(f"  n={n} k={k} e={e}: t-degrees present = {sorted(mons)}")
