"""Band certification.  At each band vertex vbar_j (j in [c0+2,e]) of split records, which of the six M0 stacks are
CERTIFIED zero by known M0 tiles, i.e. by tiles of W (M0 tiles J_o removes) lying in J_o's hexagons H^o_{j'} for j' <= j,
or in A's hexagons at steps <= p+(j-c0)?  An edge E at vbar is certified 0 when the two Lambda_k-triangles of Delta at E
have lifts that are both known and share an edge (equal heights).  Also: rhombi over the six edges in M0 (fiber lemma
hypothesis).  If the certified-zero set contains three cyclically consecutive edges, closure forces the degenerate
hexagon.  usage: python3 bandcert.py SIZE"""
import sys, os, pickle, math, cmath
from collections import Counter
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/vstar.py').read().split('st = Counter()')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
st = Counter(); fails = []
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
    rh0 = Counter(zpk(u, k) for u in M0 if u[0] == 'rh')
    sq0 = Counter(zpk(u, k) for u in M0 if u[0] == 'sq')
    for j in range(c0 + 2, e + 1):
        known = set()
        for jj in range(0, j + 1): known |= (Jo[jj][1] & M0)
        for kk in range(0, p + (j - c0) + 1):
            if kk < len(A): known |= (A[kk][1] & M0)
        overK = {zpk(u, k): u for u in known if zisL(u, k)}
        Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
        v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
        vb = zfree(v, k)
        cert = []; nrh = 0
        for i, d in enumerate(DIRS):
            nb = (vb[0] + d[0], vb[1] + d[1]); E = frozenset((vb, nb)); nrh += rh0[E]
            d2 = DIRS[(i + 1) % 6]; d0 = DIRS[(i - 1) % 6]
            t1 = frozenset((vb, nb, (vb[0] + d2[0], vb[1] + d2[1]))); t2 = frozenset((vb, nb, (vb[0] + d0[0], vb[1] + d0[1])))
            ok = False
            if t1 in overK and t2 in overK:
                ok = zfix(min(overK[t1][1]), k) == zfix(min(overK[t2][1]), k)
            # known square over E certifies NONzero
            ksq = sum(1 for u in known if u[0] == 'sq' and zpk(u, k) == E)
            cert.append('0' if ok else ('+' if ksq else '?'))
        c = ''.join(cert)
        three = any(all(c[(i + r) % 6] == '0' for r in range(3)) for i in range(6))
        st[(T, 'e' if j == e else 'sh', 'cert', c, 'three-consec-0', three, 'rh', nrh)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
