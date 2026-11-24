"""Reach of the band invariant for (A3).  For every A-shaped wake tile u read by A at step kk >= p+2 (split records):
is kk <= p+m (shadow) or after; is u a Lambda_k-triangle with a band vertex vbar_j among its projected vertices such
that proj(u) = t_{01}(vbar_j) (B) / t_{45}(vbar_j) (C), so that Theorem band gives u - alpha in M0?  And is u-alpha in M0.
usage: python3 banda3.py SIZE"""
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
    T = 'BC'[L]; Tp = OTH[T]; k = ZKK[Tp]
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1; m = e - c0
    al = (0, 0, 0, 1) if T == 'B' else (0, -1, 0, 0)
    bvs = set()
    for j in range(c0 + 2, e + 1):
        Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
        v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1))); bvs.add(zfree(v, k))
    want = (0, 1) if T == 'B' else (4, 5)
    for kk in range(p + 2, len(A)):
        for u in A[kk][1]:
            if u not in W or u[0] == 'rh': continue
            ua = (u[0], frozenset(tuple(a - b for a, b in zip(Q, al)) for Q in u[1]))
            # A-shaped: -alpha is an edge vector of u
            if not any(tuple(a - b for a, b in zip(Q2, Q1)) == tuple(-x for x in al) for Q1 in u[1] for Q2 in u[1]): continue
            cov = False
            if zisL(u, k):
                pu = zpk(u, k)
                for x in bvs:
                    if pu == tri_at(x, want[0], 1): cov = True
            st[(T, 'shadow' if kk <= p + m else 'after', u[0], 'covered by band', cov, 'u-alpha in M0', ua in M0)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
