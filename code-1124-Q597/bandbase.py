"""Base of the band induction: at j = c0+2, the edge i_T (B: 3, C: 2) of vbar_j.  Its two Lambda_k-lifts in A's contact
hexagon H_p: shape relative to o_j, whether they are A-made before p (they must not be: then they are M0 tiles by
A's no-revisit), and the type of A's step p.   usage: python3 bandbase.py SIZE"""
import sys
from collections import Counter
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/bandrule.py').read().split('st = Counter(); shapes')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTH[T]; k = ZKK[Tp]; TA = typed(A); TO = typed(Jo)
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    creA = set()
    for s in A[:p]: creA |= s[2]
    j = c0 + 2
    Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
    v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1))); vb = zfree(v, k)
    iE = SI[T][1]
    over = {zpk(u, k): u for u in A[p][1] if zisL(u, k)}
    a2, b2 = tri_at(vb, iE, 1), tri_at(vb, iE, -1)
    ua, ub = over.get(a2), over.get(b2)
    shp = tuple((x[0], tuple(sorted(tuple(c - d for c, d in zip(Q, o)) for Q in x[1]))) if x else None for x in (ua, ub))
    st[(T, 'A type p', TA[p], 'shape rel o_j', shp, 'A-made', (ua in creA, ub in creA), 'in M0', (ua in M0, ub in M0), 'in W', (ua in W, ub in W))] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
