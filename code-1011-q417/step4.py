"""Q417 STEP 4: the TARGET, from Warnaar's own formula -- INDEPENDENT of any Psi.

The refutation must NOT lean on `c_alpha(1) = A_alpha - B_alpha`.  That identity is my
own Q410 Lemma t1, and it is exactly a statement about a weight being 1 at t=1 -- using
it here would make the verification circular (memory: a-derived-statistic-can-be-blind).
So: compute c_alpha(t) from warn.py (Warnaar 2511.17034 Sect 6, G * prod x), read off the
coefficient of x^alpha, and report c_alpha(1) from THAT.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
from core import Engine
from warn import G_times_prodx
from korff import t
import sympy as sp

for (k, ell) in [(1,0),(1,1),(1,2),(2,1)]:
    M = k+ell+1; n = 2*M
    try:
        E = Engine(n)
        F = G_times_prodx(E, k, ell)
        P = sp.Poly(F, *E.x)
        alpha = (1,)*n
        ca = sp.Integer(0)
        for mon, co in P.terms():
            if tuple(mon) == alpha: ca = sp.expand(co)
        pred = sp.expand(sp.catalan(M) - (2*M-1)*t**(M-1) + t**(M+1))
        print(f"k={k} ell={ell} M={M} n={n}:")
        print(f"   WARNAAR ENGINE  c_alpha(t) = {ca}")
        print(f"   c_alpha(1) = {ca.subs(t,1)}        factored: {sp.factor(ca)}")
        print(f"   brief's closed form C_M-(2M-1)t^(M-1)+t^(M+1) = {pred}")
        print(f"   agree with closed form? {sp.expand(ca-pred)==0}  (diff {sp.expand(ca-pred)})")
    except Exception as ex:
        print(f"k={k} ell={ell}: engine failed/too big: {type(ex).__name__} {ex}")
