"""Q410 verification suite."""
import sys, itertools
sys.path.insert(0,'.')
from korff import *
from walks import AB_from_walks
from core import Engine
from warn import G_times_prodx
import sympy as sp
t = sp.Symbol('t')

def ordt1(f):
    f = sp.expand(f)
    if f == 0: return None
    k = 0
    while True:
        if sp.simplify(sp.diff(f, t, k).subs(t, 1)) != 0: return k
        k += 1
        if k > 40: return ">40"

print("=== (1) nu(T) = sum_a |J_a(T)| for alpha=(1^n): must equal n for every T ===")
for (k, ell, n) in [(1,0,4),(1,1,6),(1,2,8),(2,0,4),(2,1,6),(3,0,6)]:
    h = 2*k; w = 2*ell+2
    nus = set()
    for P in paths(h, w, tuple([1]*n)):
        nu = sum(len(J_cyclic(tuple(P[a+1][i]-P[a][i] for i in range(h))))
                 for a in range(len(P)-1))
        nus.add(nu)
    print(f"   k={k} ell={ell} n={n} (h={h},w={w}): #T={len(paths(h,w,tuple([1]*n)))}, "
          f"nu(T) values {sorted(nus)}  (predicted {{{n}}})")

print("\n=== (2) ord_{t=1} c_alpha for the witness family, M=ell+2, k=1, n=N=2M ===")
print("    c_alpha(t) = C_M - (2M-1)t^{M-1} + t^{M+1}   (Q405 thm:main)")
for M in range(2, 10):
    CM = sp.binomial(2*M, M) - sp.binomial(2*M, M-1)
    ca = CM - (2*M-1)*t**(M-1) + t**(M+1)
    print(f"   M={M}: c_alpha = {sp.expand(ca)} ; c(1)={ca.subs(t,1)} ; "
          f"ord_{{t=1}} = {ordt1(ca)} ; required ord >= n = {2*M}")

print("\n=== (3) ground truth re-derived independently (Warnaar determinant engine) ===")
for (k, ell, n) in [(1,0,4),(1,1,6)]:
    E = Engine(n); F = G_times_prodx(E, k, ell)
    ca = sp.Poly(F, *E.x).coeff_monomial(sp.prod(E.x))
    M = k+ell+1; CM = sp.binomial(2*M,M)-sp.binomial(2*M,M-1)
    pred = CM - (2*M-1)*t**(M-1) + t**(M+1)
    print(f"   k={k} ell={ell} n={n}: c_alpha={sp.expand(ca)}  matches Q405 prediction: "
          f"{sp.expand(ca-pred)==0}  ord_{{t=1}}={ordt1(ca)}")

print("\n=== (4) WHERE THE FAILURE LIVES: ordinary psi vs Korff's cyclic Psi, alpha=(1^4) ===")
print("    Warnaar Conj 1.5 at ell=0: (x..)^k G = sum_{la even, la_1<=2k} P_la(x;t).")
# ordinary HL: SSYT of shape (2,2), content (1,1,1,1), Macdonald psi via the OPEN index set
def psi_ordinary(chain):
    """chain of partitions (as conjugate-length tuples padded); Macdonald (5.11')."""
    out = sp.Integer(1)
    for a in range(len(chain)-1):
        mu, la = chain[a], chain[a+1]
        L = max(len(mu), len(la)) + 2
        mp = [sum(1 for p in mu if p >= j) for j in range(1, L+1)]
        lp = [sum(1 for p in la if p >= j) for j in range(1, L+1)]
        th = [lp[i]-mp[i] for i in range(L)]
        for i in range(L-1):
            if th[i] == 0 and th[i+1] == 1:
                out *= (1 - t**(mp[i]-mp[i+1]))
    return sp.expand(out)
chains = [[(),(1,),(2,),(2,1),(2,2)], [(),(1,),(1,1),(2,1),(2,2)]]
tot = sp.Integer(0)
for ch in chains:
    p = psi_ordinary(ch)
    nz = sum(1 for a in range(len(ch)-1)
             if psi_ordinary([ch[a],ch[a+1]]) != 1)
    print(f"   ordinary SSYT chain {ch}: psi_T = {sp.factor(p)}  "
          f"(#vanishing factors = ord_{{t=1}} = {ordt1(p)} out of 4 strips)")
    tot += p
print(f"   SUM over the 2 ordinary SSYT of shape (2,2) = {sp.factor(tot)} = {sp.expand(tot)}")
print(f"   target c_alpha(t) = 2 - 3t + t^3 : MATCH = {sp.expand(tot-(2-3*t+t**3))==0}")
print("   Korff's cyclic weight on the 4 CYLINDRIC tableaux of content (1^4), h=w=2:")
for P in paths(2,2,(1,1,1,1)):
    nu = sum(len(J_cyclic(tuple(P[a+1][i]-P[a][i] for i in range(2)))) for a in range(4))
    print(f"      {P[-1]}  Psi_T = {sp.factor(weight_path(P,2))}  nu={nu}  ord_{{t=1}}={ordt1(weight_path(P,2))}")
print("   => cyclic index set forces a vanishing factor at EVERY one of the 4 strips (ord 4);")
print("      the open index set forces one at exactly 2 of them (ord 2).  2 != 4.")
