"""(i) Locus as a finite configuration statement: for each A-shaped wake tile u read by A at kk = p+i (2 <= i <= m),
its shape relative to o = min R^o_{c0+i} and the type of J_o's step c0+i (= the shadow hexagon H_{p+i} = H^o_{c0+i}?).
(ii) (H4) stack chains: at band vertices, the squares over E_1 and over E_4 have Lambda_{1-k}-edges exactly the k1
unit segments of the line from h(G1) to h(G2).   usage: python3 banda3c.py SIZE"""
import sys
from collections import Counter, defaultdict
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/bandrule.py').read().split('st = Counter(); shapes')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
st = Counter(); shp = defaultdict(Counter)
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
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1; m = e - c0
    al = (0, 0, 0, 1) if T == 'B' else (0, -1, 0, 0)
    for i in range(2, m + 1):
        st[(T, 'H_{p+i} == H^o_{c0+i}', A[p + i][1] == Jo[c0 + i][1])] += 1
        o = min(next(u for u in Jo[c0 + i][1] if u[0] == 'rh')[1])
        for u in A[p + i][1]:
            if u not in W or u[0] == 'rh': continue
            if not any(tuple(a - b for a, b in zip(Q2, Q1)) == tuple(-x for x in al) for Q1 in u[1] for Q2 in u[1]): continue
            shp[(T, TO[c0 + i])][(u[0], tuple(sorted(tuple(a - b for a, b in zip(Q, o)) for Q in u[1])))] += 1
    # (H4)
    SQ = defaultdict(set)
    for u in M0:
        if u[0] == 'sq': SQ[zpk(u, k)].add(frozenset(zfix(Q, k) for Q in u[1]))
    lift = {zpk(u, k): zfix(min(u[1]), k) for u in M0 if zisL(u, k)}
    for j in range(c0 + 2, e + 1):
        o = min(next(u for u in Jo[j][1] if u[0] == 'rh')[1])
        v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1))); x = zfree(v, k)
        hG1, hG2 = lift[tri_at(x, 2, 1)], lift[tri_at(x, 5, 1)]
        n = (hG1[0] - hG2[0], hG1[1] - hG2[1]); kk1 = max(abs(n[0]), abs(n[1]))
        segs = set()
        if kk1:
            un = (n[0] // kk1, n[1] // kk1)
            segs = {frozenset(((hG2[0] + r * un[0], hG2[1] + r * un[1]), (hG2[0] + (r + 1) * un[0], hG2[1] + (r + 1) * un[1]))) for r in range(kk1)}
        E1 = frozenset((x, (x[0], x[1] + 1))); E4 = frozenset((x, (x[0], x[1] - 1)))
        st[(T, '(H4) chains', SQ[E1] == segs and SQ[E4] == segs)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
for kk in sorted(shp, key=str):
    print('LOCUS', kk, dict(shp[kk]))
