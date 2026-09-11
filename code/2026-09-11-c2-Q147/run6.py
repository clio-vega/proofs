"""Split the commutator by bead-number sector and solve each locus separately."""
import sympy as sp
from engine import partitions, commutator, maya, unmaya, apply_R

def sector_eqs(e, f, maxn, shared=True, NM=10):
    if shared:
        W = [sp.Integer(1)] + [sp.Symbol('w%d'%i) for i in range(1,NM)]
        tab = {e: W, f: W}
    else:
        tab = {e: [sp.Integer(1)] + [sp.Symbol('w%d'%i) for i in range(1,e)],
               f: [sp.Integer(1)] + [sp.Symbol('W%d'%i) for i in range(1,f)]}
    w = lambda g,N: tab[g][N]
    one, two = set(), set()
    for n in range(maxn+1):
        for lam in partitions(n):
            L = n + e + f + 4
            M = maya(lam, L)
            for mu, val in commutator(lam, e, f, w).items():
                Mp = maya(mu, L)
                d = len(M ^ Mp)
                (one if d == 2 else two).add(sp.expand(val))
    return one, two

for (e,f,mx) in [(2,3,11),(2,4,11),(3,4,12),(2,5,12),(3,5,13)]:
    one, two = sector_eqs(e,f,mx)
    syms = sorted({s for q in one|two for s in q.free_symbols}, key=str)
    s1 = sp.solve(list(one), syms, dict=True)
    s2 = sp.solve(list(two), syms, dict=True)
    print("e,f=(%d,%d) shared: 1-bead #%d -> %s" % (e,f,len(one),s1))
    print("            2-bead #%d -> %s" % (len(two), s2))
