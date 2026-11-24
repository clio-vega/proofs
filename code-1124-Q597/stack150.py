"""Provenance of the squares stacked over the edges of vbar (J_o's board before step e, split records): for each edge
angle (from tau3), the list of square provenances (creJ@d, W(removed at e+d), X) and whether the square lies in H*.
usage: python3 stack150.py SIZE"""
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
    T = 'BC'[L]; Tp = OTH[T]; k = ZKK[Tp]
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    rem = {}; creJ = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
        for u in s[2]: creJ.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    bd = set(M0)
    for s in Jo[:e]: bd = (bd - s[1]) | s[2]
    Hs = Jo[e][1]; Rs = next(u for u in Hs if u[0] == 'rh'); o = min(Rs[1])
    v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
    vb = zfree(v, k)
    t3 = next(u for u in Hs if zisL(u, k) and v in u[1]); c = [x for x in zpk(t3, k) if x != vb]
    a3 = ang(vb, ((c[0][0] + c[1][0]) / 2, (c[0][1] + c[1][1]) / 2))
    out = []
    for d in [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]:
        nb = (vb[0] + d[0], vb[1] + d[1]); edge = frozenset((vb, nb))
        sqs = [u for u in bd if u[0] == 'sq' and zpk(u, k) == edge]
        if not sqs: continue
        pv = tuple(sorted(('H*' if u in Hs else 'creJ@%+d' % (creJ[u] - e) if u in creJ else 'W@%+d' % (rem[u] - e) if u in W else 'X') for u in sqs))
        out.append((round((ang(vb, nb) - a3) % 360), pv))
    st[('T', T, tuple(sorted(out)))] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
