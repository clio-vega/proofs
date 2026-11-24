"""Locality of the band induction.  At band vertex j, known = M0 tiles of J_o's hexagons H^o_{j+a} for a in WIN,
plus A's contact hexagon H_p if USEHP, plus (IH) the M0 lifts at vbar_{j-1} only.  Also records which certified zero
edges come from which source.  usage: python3 bandloc.py SIZE WIN USEHP   (WIN like -1,0)"""
import sys, itertools
from collections import Counter
SIZE = int(sys.argv[1]); WIN = [int(x) for x in sys.argv[2].split(',')]; USEHP = int(sys.argv[3])
exec(open('/home/clio/projects/proofs/code-1124-Q597/vstar.py').read().split('st = Counter()')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
SOLS = [K for K in itertools.product(range(7), repeat=6) if K[0]-K[3] == K[4]-K[1] == K[2]-K[5]]
def forced(c, idxs):
    ok = [K for K in SOLS if not any(c[i] == '0' and K[i] for i in range(6)) and not any(c[i] == '+' and not K[i] for i in range(6))]
    return all(all(K[i] == 0 for i in idxs) for K in ok)
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
    over0 = {zpk(u, k): u for u in M0 if zisL(u, k)}
    prevlifts = {}; prev = None; allok = True
    for j in range(c0 + 2, e + 1):
        known = set()
        for a in WIN:
            if 0 <= j + a <= e: known |= (Jo[j + a][1] & M0)
        if USEHP == 1 or (USEHP == 2 and j == c0 + 2): known |= (A[p][1] & M0)
        overK = {zpk(u, k): u for u in known if zisL(u, k)}
        src = {tau: 'J' for tau in overK}
        for tau, u in prevlifts.items():
            if tau not in overK: overK[tau] = u; src[tau] = 'IH'
        sqK = Counter(zpk(u, k) for u in known if u[0] == 'sq')
        Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
        v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
        vb = zfree(v, k)
        c = ''; how = []
        for i, d in enumerate(DIRS):
            nb = (vb[0] + d[0], vb[1] + d[1]); E = frozenset((vb, nb))
            t1 = frozenset((vb, nb, (vb[0] + DIRS[(i+1) % 6][0], vb[1] + DIRS[(i+1) % 6][1])))
            t2 = frozenset((vb, nb, (vb[0] + DIRS[(i-1) % 6][0], vb[1] + DIRS[(i-1) % 6][1])))
            if t1 in overK and t2 in overK and zfix(min(overK[t1][1]), k) == zfix(min(overK[t2][1]), k):
                c += '0'; how.append((i, src[t1], src[t2]))
            elif sqK[E] or (t1 in overK and t2 in overK): c += '+'
            else: c += '?'
        ok = forced(c, (0, 2, 3, 5)); allok &= ok
        step = None if prev is None else (vb[0] - prev[0], vb[1] - prev[1])
        st[(T, 'j-c0', min(j - c0, 3), TO[j], 'step', step, 'cert', c, 'ok', ok, 'zeros-from', tuple(how))] += 1
        prevlifts = {tau: u for tau, u in over0.items() if vb in tau}; prev = vb
    st[(T, 'RECORD', allok)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
