"""v_S(Atilde_{I,I}) at w_m = eps, measured per I, two mechanisms:
   (1) the engine's Atilde (symmetrisation / E^X quotient)
   (2) diag_rootproduct (independent: product over positive roots of C_N)
DOMAIN: N=2,3; X=C; m in [N]; eps in {+1,-1}; q0=3/7; other w_j random rationals.
"""
import sys, os, random
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code-1017c2-coboundary'))
import q484_engine as E
t = sp.Symbol('t')

def ordt(expr):
    expr = sp.cancel(sp.together(expr))
    if expr == 0:
        return None
    num, den = sp.fraction(expr)
    pn, pd = sp.Poly(sp.expand(num), t), sp.Poly(sp.expand(den), t)
    def lo(p):
        c = list(reversed(p.all_coeffs()))
        for i, x in enumerate(c):
            if sp.simplify(x) != 0:
                return i
        return None
    a, b = lo(pn), lo(pd)
    return None if a is None else a - b

for N in (2, 3):
  for m in range(1, N+1):
    for eps in (1, -1):
        random.seed(5)
        w = []
        for j in range(1, N+1):
            w.append(sp.Integer(eps)+t if j == m else sp.Rational(random.randint(3,40), random.randint(1,7)))
        q0 = sp.Rational(3,7)
        tot_e = tot_r = 0
        rows = []
        for I in E.subsets(N):
            ae = E.Atilde(N, w, I, I, X='C').subs({E.q: q0})
            ar = E.diag_rootproduct(N, w, I).subs({E.q: q0})
            oe, orr = ordt(ae), ordt(ar)
            tot_e += oe or 0; tot_r += orr or 0
            rows.append((I, oe, orr, m in I))
        print(f"N={N} m={m} eps={eps}: sum_engine={tot_e} sum_rootprod={tot_r}  #{{I:m in I}}={sum(1 for I in E.subsets(N) if m in I)}")
        for I, oe, orr, inI in rows:
            flag = '' if oe == (1 if inI else 0) else '   <-- not [m in I]'
            print(f"    I={str(I):12s} v_engine={oe} v_rootprod={orr}  m in I={inI}{flag}")
