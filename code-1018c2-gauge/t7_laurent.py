"""Is Atil_{M,I} itself a LAURENT POLYNOMIAL?  (i.e. denominator a monomial)

If yes then every height-one prime of R=C(q)[w^pm] at which Atil has a pole is
EMPTY, and the divisor-by-divisor argument needs only the divisors of det Phi.
Also: the second orientation of the (Q2) family, {1 - q w_1/w_2 = 0} (e>c),
which the first-pass proof of v(F)>=1 does NOT cover.
"""
import sympy as sp, itertools, sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1018c2-gauge')
from q484_engine import q, wvars, Atilde

def subsets(N):
    return [tuple(s) for k in range(N+1) for s in itertools.combinations(range(1,N+1),k)]

def is_monomial(expr, gens):
    """True iff expr is c * monomial in gens."""
    e = sp.expand(expr)
    try:
        p = sp.Poly(e, *gens)
    except sp.PolynomialError:
        return False
    return len(p.terms()) == 1

for N in (1,2,3):
    w = wvars(N)
    gens = [q] + list(w)
    ok = n = 0
    bad = []
    for M in subsets(N):
        for I in subsets(N):
            a = Atilde(N, w, M, I, X='C')
            if a == 0: continue
            num, den = sp.fraction(sp.cancel(sp.together(a)))
            n += 1
            if is_monomial(den, gens): ok += 1
            else: bad.append((M, I, sp.factor(den)))
    print(f"N={N}: Atilde Laurent (denominator a monomial): {ok}/{n}")
    for b in bad[:4]:
        print(f"    NOT Laurent: M={b[0]} I={b[1]}  den={b[2]}")
