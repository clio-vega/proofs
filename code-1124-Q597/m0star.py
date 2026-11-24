"""Fiber star at vbar in M0 (not J_o's board): the 6 Delta-triangles around vbar with their M0 heights (rel. o's fixed part)
and provenance (W removed by J_o at e+d / X), and the stacks over the 6 edges with each square's provenance, plus rhombi.
Split records, T in {B, C}.   usage: python3 m0star.py SIZE"""
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
    rem = {}; remA = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    for kk, s in enumerate(A):
        for u in s[1]: remA.setdefault(u, kk)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    Hs = Jo[e][1]; Rs = next(u for u in Hs if u[0] == 'rh'); o = min(Rs[1])
    v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
    vb = zfree(v, k); h0 = zfix(v, k)
    t3 = next(u for u in Hs if zisL(u, k) and v in u[1]); c = [x for x in zpk(t3, k) if x != vb]
    a3 = ang(vb, ((c[0][0] + c[1][0]) / 2, (c[0][1] + c[1][1]) / 2))
    def prov(u):
        s = 'W@e%+d' % (rem[u] - e) if u in rem else 'X'
        if u in remA: s += '/A@p%+d' % (remA[u] - p)
        return s
    over = {zpk(u, k): u for u in M0 if zisL(u, k)}
    tris = []
    for tau, u in over.items():
        if vb in tau:
            cc = [x for x in tau if x != vb]; mid = ((cc[0][0] + cc[1][0]) / 2, (cc[0][1] + cc[1][1]) / 2)
            tris.append((round((ang(vb, mid) - a3) % 360), tuple(int(a - b) for a, b in zip(zfix(min(u[1]), k), h0)), prov(u)))
    stacks = []
    for d in [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]:
        nb = (vb[0] + d[0], vb[1] + d[1]); edge = frozenset((vb, nb))
        its = [u for u in M0 if u[0] in ('sq', 'rh') and zpk(u, k) == edge]
        if its: stacks.append((round((ang(vb, nb) - a3) % 360), tuple(sorted(u[0] + ':' + prov(u) for u in its))))
    nfib = sum(1 for u in M0 if u[0] == 'tri' and len({zfree(Q, k) for Q in u[1]}) == 1 and zfree(next(iter(u[1])), k) == vb)
    st[('T', T, 'tris', tuple(sorted(tris)), 'stacks', tuple(sorted(stacks)), '#fiber tris', nfib)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
