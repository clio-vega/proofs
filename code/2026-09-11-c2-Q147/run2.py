import sympy as sp
from engine import partitions, commutator

NMAX = 10
W = [sp.Integer(1)] + [sp.Symbol('w%d' % i) for i in range(1, NMAX)]
def wshared(e, N):  return W[N]

def collect(e, f, maxn):
    eqs = set()
    for n in range(0, maxn + 1):
        for lam in partitions(n):
            for mu, val in commutator(lam, e, f, wshared).items():
                eqs.add(sp.expand(val))
    return eqs

for (e, f, mx) in [(2,4,11),(3,4,12),(2,5,12),(3,5,13),(4,5,13),(1,5,10),(2,6,13)]:
    eqs = collect(e, f, mx)
    syms = [W[i] for i in range(1, NMAX)]
    sol = sp.solve(list(eqs), syms, dict=True)
    used = sorted({s for q in eqs for s in q.free_symbols}, key=str)
    print("e,f=(%d,%d) maxn=%d  #eqs=%d  symbols seen=%s" % (e,f,mx,len(eqs), used))
    print("    SOLUTIONS:", sol)
