import sympy as sp
from engine import partitions, commutator

def mk(e):  # w^{(e)}_0 = 1, rest free
    return [sp.Integer(1)] + [sp.Symbol('v%d_%d' % (e, i)) for i in range(1, e)]

def run(e, f, maxn):
    we, wf = mk(e), mk(f)
    tab = {e: we, f: wf}
    w = lambda g, N: tab[g][N]
    eqs = set()
    for n in range(0, maxn + 1):
        for lam in partitions(n):
            for mu, val in commutator(lam, e, f, w).items():
                eqs.add(sp.expand(val))
    syms = [s for s in we[1:] + wf[1:]]
    sol = sp.solve(list(eqs), syms, dict=True)
    print("e,f=(%d,%d) maxn=%d #eqs=%d syms=%s" % (e,f,maxn,len(eqs), syms))
    print("   SOLUTIONS:", sol)
    return eqs

for (e,f,mx) in [(1,3,9),(2,3,11),(2,4,11),(3,4,12),(2,5,12),(3,5,13),(4,5,13),(3,6,14)]:
    run(e,f,mx)
