import sympy as sp
from ribbon import *
from itertools import product
t = sp.Symbol('t')

def compositions(n):
    if n==0: yield (); return
    for k in range(1,n+1):
        for rest in compositions(n-k): yield (k,)+rest

def apply_R(vec, e, N):
    """vec: dict lam->coeff.  Returns R_e(t) vec."""
    out={}
    for mu,c in vec.items():
        for lam,ht in add_ribbons(mu,e,N+3):
            out[lam]=sp.expand(out.get(lam,0)+c*t**ht)
    return out

def g_comp(beta, N):
    """R_{beta_l} ... R_{beta_1} . 1   (beta_1 added FIRST)"""
    v={():1}
    for e in beta: v=apply_R(v,e,N)
    return v

print("=== T3a: do R_a(t), R_b(t) commute?  (sorting condition) ===")
for (a,b) in [(1,2),(1,3),(2,3),(2,4),(1,4)]:
    n=a+b
    v1=g_comp((a,b),n); v2=g_comp((b,a),n)
    keys=set(v1)|set(v2)
    diff={k:sp.factor(sp.expand(v1.get(k,0)-v2.get(k,0))) for k in keys}
    diff={k:v for k,v in diff.items() if v!=0}
    print(f"  R_{b}R_{a}1 - R_{a}R_{b}1  (n={n}):", diff if diff else "ZERO")

print()
print("=== T3b: for which t does the sorting condition hold? ===")
allroots=set()
for (a,b) in [(1,2),(1,3),(2,3),(2,4),(3,4),(1,4),(2,5)]:
    n=a+b
    v1=g_comp((a,b),n); v2=g_comp((b,a),n)
    for k in set(v1)|set(v2):
        d=sp.expand(v1.get(k,0)-v2.get(k,0))
        if d!=0:
            rs=sp.solve(sp.Eq(d,0),t)
            allroots.add((a,b,k,tuple(sorted(map(str,rs)))))
sols=None
for (a,b,k,rs) in allroots:
    S=set(rs)
    sols = S if sols is None else (sols & S)
print("  common roots over all tested (a,b) and all Schur coefficients:", sols)
