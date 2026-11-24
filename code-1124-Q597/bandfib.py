"""Band fibers in M0.  For every shadow step j in [c0+2, e] of split records: vertex v_j (same rule as at the split),
in M0: number of fiber triangles over vbar_j (other-lattice triangles whose free part is vbar_j), the stack sizes on the 6
edges at vbar_j (angle from tau3_j = the triangle of H^o_j at v_j of the projected lattice, if any), and whether
the 150-edge stack is empty.  usage: python3 bandfib.py SIZE"""
import sys, os, pickle, math, cmath
from collections import Counter
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/vstar.py').read().split('st = Counter()')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTH[T]; k = ZKK[Tp]; TO = typed(Jo)
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    over = {zpk(u, k): u for u in M0 if zisL(u, k)}
    sq = Counter(zpk(u, k) for u in M0 if u[0] == 'sq')
    fib = Counter(zfree(next(iter(u[1])), k) for u in M0 if u[0] == 'tri' and len({zfree(Q, k) for Q in u[1]}) == 1)
    for j in range(c0 + 2, e + 1):
        Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
        v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
        vb = zfree(v, k)
        t3 = [u for u in Jo[j][1] if zisL(u, k) and v in u[1]]
        ks = []
        for d in [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]:
            ks.append(sq[frozenset((vb, (vb[0] + d[0], vb[1] + d[1])))])
        st[('T', T, 'step', 'e' if j == e else 'sh', TO[j], 'tau3 present', len(t3), '#fiber', fib[vb], 'stacks', tuple(ks))] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
