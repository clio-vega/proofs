"""Induction along the band.  For j = c0+2..e: known Lambda_k-lifts = lifts in known M0 tiles (J_o hexagons <= e, A's
hexagons <= p+m) PLUS (induction hypothesis) the M0 lifts of all triangles of Delta at vbar_{j'} for band j' < j.
Certify stacks at vbar_j (0 / + / ?), solve closure, and report whether k0=k2=k3=k5=0 is forced (degenerate on axis
{1,4}).  Also the band step direction vbar_j - vbar_{j-1}.  usage: python3 bandind.py SIZE"""
import sys, itertools
from collections import Counter
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/vstar.py').read().split('st = Counter()')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
SOLS = [K for K in itertools.product(range(7), repeat=6) if K[0]-K[3] == K[4]-K[1] == K[2]-K[5]]
def forced(c, idxs):
    ok = [K for K in SOLS if not any(c[i] == '0' and K[i] for i in range(6)) and not any(c[i] == '+' and not K[i] for i in range(6))]
    return all(all(K[i] == 0 for i in idxs) for K in ok)
st = Counter(); bad = []
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
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    over0 = {zpk(u, k): u for u in M0 if zisL(u, k)}
    known = set()
    for jj in range(c0, e + 1): known |= (Jo[jj][1] & M0)  # VARIANT: J_o from contact
    known |= (A[p][1] & M0)  # VARIANT: A contact hexagon only
    overK = {zpk(u, k): u for u in known if zisL(u, k)}
    sqK = Counter(zpk(u, k) for u in known if u[0] == 'sq')
    prev = None; allok = True
    for j in range(c0 + 2, e + 1):
        Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
        v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
        vb = zfree(v, k)
        step = None if prev is None else (vb[0] - prev[0], vb[1] - prev[1])
        c = ''
        for i, d in enumerate(DIRS):
            nb = (vb[0] + d[0], vb[1] + d[1]); E = frozenset((vb, nb))
            t1 = frozenset((vb, nb, (vb[0] + DIRS[(i+1) % 6][0], vb[1] + DIRS[(i+1) % 6][1])))
            t2 = frozenset((vb, nb, (vb[0] + DIRS[(i-1) % 6][0], vb[1] + DIRS[(i-1) % 6][1])))
            if t1 in overK and t2 in overK and zfix(min(overK[t1][1]), k) == zfix(min(overK[t2][1]), k): c += '0'
            elif sqK[E]: c += '+'
            elif t1 in overK and t2 in overK: c += '+'   # both lifts known, different heights
            else: c += '?'
        ok = forced(c, (0, 2, 3, 5))
        allok &= ok
        st[(T, 'e' if j == e else 'sh', 'step', step, 'cert', c, 'degenerate forced', ok)] += 1
        # induction hypothesis: add all lifts at vbar_j
        for tau, u in over0.items():
            if vb in tau: overK[tau] = u
        prev = vb
    st[(T, 'RECORD all band vertices forced', allok)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
