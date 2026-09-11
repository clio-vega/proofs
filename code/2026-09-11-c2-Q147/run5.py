import sympy as sp
from engine import partitions, commutator

# --- (i) does case (B) of the one-bead sector alone already force the anchor?
def caseB_eqs(e, f):
    W  = [sp.Integer(1)] + [sp.Symbol('w%d'%i) for i in range(1,e)]
    WB = [sp.Integer(1)] + [sp.Symbol('W%d'%i) for i in range(1,f)]
    eqs = []
    for u in range(e):
        for v in range(f-e):
            for z in range(e):
                eqs.append(sp.expand(W[z]*WB[u+v+1] + W[u]*WB[v+z]))
    syms = [s for s in W[1:]+WB[1:]]
    return sp.solve(eqs, syms, dict=True), syms

for (e,f) in [(1,2),(1,3),(1,5),(2,3),(3,4),(4,5),(2,5),(3,7),(5,6),(4,9),(6,7)]:
    sol, syms = caseB_eqs(e,f)
    ok = (len(sol)==1 and all(sol[0][s]==(-1)**int(str(s)[1:]) for s in syms))
    print("case(B) alone, e,f=(%d,%d): unique anchor? %s   %s" % (e,f,ok, sol if not ok else ''))

# --- (ii) drop the normalisation w_0 = 1 entirely (shared weights)
print()
NM=8
V=[sp.Symbol('w%d'%i) for i in range(NM)]
def w(e,N): return V[N]
for (e,f,mx) in [(1,2,7),(1,3,8),(2,3,11),(2,4,11)]:
    eqs=set()
    for n in range(mx+1):
        for lam in partitions(n):
            for mu,val in commutator(lam,e,f,w).items(): eqs.add(sp.expand(val))
    used=sorted({s for q in eqs for s in q.free_symbols},key=str)
    sol=sp.solve(list(eqs), used, dict=True)
    print("unnormalised e,f=(%d,%d): syms=%s  SOL=%s" % (e,f,used,sol))
