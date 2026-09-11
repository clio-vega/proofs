import sympy as sp
from engine import partitions, commutator

NMAX = 8
W = [sp.Integer(1)] + [sp.Symbol('w%d' % i) for i in range(1, NMAX)]
def wshared(e, N):  return W[N]

def collect_eqs(e, f, maxn, w=wshared):
    eqs = set()
    for n in range(0, maxn + 1):
        for lam in partitions(n):
            for mu, val in commutator(lam, e, f, w).items():
                eqs.add(sp.expand(val))
    return eqs

for (e, f) in [(1,2),(1,3),(1,4),(2,3)]:
    eqs = collect_eqs(e, f, 7)
    print("=== e,f =", e, f, " #distinct nonzero commutator entries:", len(eqs))
    for q in sorted(eqs, key=lambda z: (sp.count_ops(z), str(z)))[:12]:
        print("   ", q)
    sol = sp.solve(list(eqs), [W[i] for i in range(1, NMAX)], dict=True)
    print("   SOLUTIONS:", sol)
