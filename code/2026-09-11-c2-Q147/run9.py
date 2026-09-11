"""The degenerate stratum: drop the normalisation w_0 = 1 (independent weights)."""
import sympy as sp
from engine import partitions, commutator

def locus(e, f, maxn):
    W  = [sp.Symbol('a%d'%i) for i in range(e)]
    WB = [sp.Symbol('b%d'%i) for i in range(f)]
    tab = {e: W, f: WB}; w = lambda g,N: tab[g][N]
    eqs = set()
    for n in range(maxn+1):
        for lam in partitions(n):
            for mu,v in commutator(lam,e,f,w).items(): eqs.add(sp.expand(v))
    syms = W + WB
    return eqs, syms, W, WB

for (e,f,mx) in [(2,3,11),(2,4,11),(3,4,12)]:
    eqs, syms, W, WB = locus(e,f,mx)
    print("=== e,f=(%d,%d)" % (e,f))
    print("   full locus:", sp.solve(list(eqs), syms, dict=True))
    # stratum a0 = 0 :
    E0 = {sp.expand(q.subs(W[0],0)) for q in eqs} - {0}
    print("   stratum a0=0 ->", sp.solve(list(E0), syms[1:], dict=True))
    E1 = {sp.expand(q.subs(WB[0],0)) for q in eqs} - {0}
    print("   stratum b0=0 ->", sp.solve(list(E1), [s for s in syms if s is not WB[0]], dict=True))
