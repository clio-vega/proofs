r"""FULL-q DIRECT CHECK of hyp:L^S, extending the measured domain to N=5.

This is INDEPENDENT of degen.py: it reads ord_t of the full-q entry at
w_m = eps + t with the other w's and q at a gated-generic rational point,
rather than taking a q=0 degeneration limit.  Two mechanisms, one claim.

DOMAIN: N = 5 (and 1..4 for regression); X in {B,C}; every m in [N];
        eps = +1, -1 (C) / +1 (B); every (A,B) with m notin A.
        One gated-generic rational (q, w) point per (N,X,m,eps).
GATE on the sample point (brief instrument fault 2): require
        ord_t diag_rootproduct(N,w,I) == [m in I]  for EVERY I,
i.e. (BD) holds at the drawn point.  Reject and redraw otherwise.
PREDICTION, written before the run:  ord_t Atilde_{A,B} == 0  for every A >= B
with m notin A;  Atilde_{A,B} == 0 identically for A not>= B (plant).
"""
import sys, os, random, argparse
import sympy as sp
sys.path.insert(0, '.')
import q484_engine as E
t = sp.Symbol('t')

def ordt(expr):
    expr = sp.cancel(sp.together(expr))
    if expr == 0:
        return None
    num, den = sp.fraction(expr)
    def lo(p):
        c = list(reversed(sp.Poly(sp.expand(p), t).all_coeffs()))
        for i, x in enumerate(c):
            if x != 0:
                return i
        return None
    a, b = lo(num), lo(den)
    return None if a is None else a - b

def generic_point(N, m, eps, q0, seed):
    random.seed(seed)
    sub = {E.q: sp.Rational(q0)}
    for _ in range(400):
        w = [sp.Integer(eps) + t if j == m else
             sp.Rational(random.randint(11, 400), random.randint(1, 23))
             for j in range(1, N + 1)]
        if all(ordt(E.diag_rootproduct(N, w, I).subs(sub)) == (1 if m in I else 0)
               for I in E.subsets(N)):
            return w
    raise RuntimeError('no generic point found in 400 draws')

ap = argparse.ArgumentParser()
ap.add_argument('--N', type=int, required=True)
ap.add_argument('--q0', default='3/11')
ap.add_argument('--seed', type=int, default=20261020)
a = ap.parse_args()
N = a.N
sub = {E.q: sp.Rational(a.q0)}
P = E.subsets(N)
good = badpos = plant_ok = plant_bad = 0
fails = []
for X in ('C', 'B'):
    for eps in ((1, -1) if X == 'C' else (1,)):
        for m in range(1, N + 1):
            w = generic_point(N, m, eps, a.q0, a.seed + 97 * m + (0 if eps == 1 else 7))
            for A in P:
                if m in A: continue
                for B in P:
                    v = E.Atilde(N, w, A, B, X=X, simplify=False).subs(sub)
                    o = ordt(v)
                    if E.geq(A, B):
                        if o == 0: good += 1
                        else:
                            badpos += 1
                            if len(fails) < 10: fails.append(('L^S FAIL', X, m, eps, A, B, o))
                    else:
                        if o is None: plant_ok += 1
                        else:
                            plant_bad += 1
                            if len(fails) < 10: fails.append(('PLANT FAIL', X, m, eps, A, B, o))
            print(f"  done X={X} eps={eps:+d} m={m}: cum good={good} badpos={badpos} "
                  f"plant_ok={plant_ok} plant_bad={plant_bad}", flush=True)
print(f"\nN={N}  q={a.q0}  seed={a.seed}")
print(f"  hyp:L^S  (ord_t == 0 required): {good} ok, {badpos} FAIL")
print(f"  plant (A not>= B, must vanish): {plant_ok} ok, {plant_bad} FAIL")
for f in fails: print("   ", f)
