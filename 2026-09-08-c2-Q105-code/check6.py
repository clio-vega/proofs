"""CHECK 6: the zero/pole decomposition.
  (z1 - b z2) Htilde_{a,b}(z1) Htilde_{c,d}(z2) = (z1 - a z2) Phi',
  (z2 - d z1) Htilde_{c,d}(z2) Htilde_{a,b}(z1) = (z2 - c z1) Phi',
with Phi' = sigma[Xz1]sigma[Xz2] f[X - (a-b)/z1 - (c-d)/z2] normal-ordered.
Four INDEPENDENT symbols a,b,c,d.  Zero of the contraction = a; pole = b.
"""
import sympy as sp
from vertex import *
a, b, c, d = sp.symbols('a b c d')

def shift2(f, cf1, cf2):
    """f[X - c1/z1 - c2/z2] -> dict (i,i') -> element."""
    acc = {(0, 0, ()): sp.Integer(1)}
    out = {}
    for mu, coef in f.items():
        acc = {(0, 0, ()): sp.Integer(1)}
        for k in mu:
            nxt = {}
            for (i, j, nu), cc in acc.items():
                for key, val in ((  (i, j, tuple(sorted(nu+(k,), reverse=True))), cc),
                                 (  (i+k, j, nu), -cc*cf1(k)),
                                 (  (i, j+k, nu), -cc*cf2(k))):
                    nxt[key] = sp.expand(nxt.get(key, 0) + val)
            acc = nxt
        for (i, j, nu), cc in acc.items():
            out.setdefault((i, j), {})
            out[(i, j)][nu] = sp.expand(out[(i, j)].get(nu, 0) + coef*cc)
    return {k: {p: v for p, v in e.items() if sp.expand(v) != 0} for k, e in out.items()}

def Phi(f, m, n, cf1, cf2):
    res = {}
    for (i, j), elt in shift2(f, cf1, cf2).items():
        pr = pmul(pmul(h(m+i), h(n+j)), elt)
        for k, v in pr.items():
            res[k] = sp.expand(res.get(k, 0) + v)
    return {k: v for k, v in res.items() if sp.expand(v) != 0}

cf1 = lambda k: a**k - b**k
cf2 = lambda k: c**k - d**k

def X(f, m, n):   # Htilde_{a,b}(z1) Htilde_{c,d}(z2), coeff of z1^m z2^n
    return Hmode(Hmode(f, n, cf2), m, cf1)
def Y(f, m, n):   # Htilde_{c,d}(z2) Htilde_{a,b}(z1)
    return Hmode(Hmode(f, m, cf1), n, cf2)

ok1 = ok2 = bad = 0
for deg in range(0, 3):
    for mu in parts(deg):
        f = {mu: sp.Integer(1)}
        for m in range(-1, 3):
            for n in range(-1, 3):
                P = Phi(f, m, n, cf1, cf2)
                Pm = Phi(f, m, n+1, cf1, cf2)   # for the z2-shifted term
                # (z1-b z2)*XX at z1^m z2^n  =  X_{m-1,n} - b X_{m,n-1}
                L1 = sub(X(f, m-1, n), smul(b, X(f, m, n-1)))
                R1 = sub(Phi(f, m-1, n, cf1, cf2), smul(a, Phi(f, m, n-1, cf1, cf2)))
                D1 = {k: sp.simplify(v) for k, v in sub(L1, R1).items()}
                D1 = {k: v for k, v in D1.items() if v != 0}
                if D1: bad += 1; print("FAIL-1", mu, m, n, D1)
                else: ok1 += 1
                # (z2-d z1)*YY  =  Y_{m,n-1} - d Y_{m-1,n}   vs  (z2-c z1)Phi
                L2 = sub(Y(f, m, n-1), smul(d, Y(f, m-1, n)))
                R2 = sub(Phi(f, m, n-1, cf1, cf2), smul(c, Phi(f, m-1, n, cf1, cf2)))
                D2 = {k: sp.simplify(v) for k, v in sub(L2, R2).items()}
                D2 = {k: v for k, v in D2.items() if v != 0}
                if D2: bad += 1; print("FAIL-2", mu, m, n, D2)
                else: ok2 += 1
print("zero/pole decomposition, 4 free symbols:", ok1, "+", ok2, "agree,", bad, "mismatch")
