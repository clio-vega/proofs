"""Height-function picture.  On J_o's board BEFORE step e (split records): for tau in {tau1 = pi1(Z1), tau2 = pi1(Z2),
tau3 = pi1(H*'s L1 triangle at v)}: the three Delta1-neighbours (unit L1 triangles sharing an edge), the height h(nbr) - o0
(Lambda0 part) of the unique L1 triangle of the board over each neighbour, and its provenance (in H*, J_o-created at step
j-e, A-created, W / X), plus whether the stack over the shared edge contains a rhombus.
usage: python3 nbrs.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1])
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0'); OTH = {'B': 'C', 'C': 'B'}
Zoff = {'C': [((0, 0, 0, -1), (0, 0, 0, 0), (0, 0, 1, -1)), ((0, 0, 0, 0), (0, 0, 1, -1), (0, 0, 1, 0))],
        'B': [((0, 1, 0, 1), (0, 2, 0, 1), (1, 1, 0, 1)), ((0, 1, 0, 1), (1, 0, 0, 1), (1, 1, 0, 1))]}
ZKK = {'C': 1, 'B': 0}   # J_o->C: Z are L1-triangles (fixed Lambda0 part); J_o->B: L0-triangles
def zfix(P, k): return P[:2] if k == 1 else P[2:]
def zfree(P, k): return P[2:] if k == 1 else P[:2]
def zisL(t, k): return t[0] == 'tri' and len({zfix(P, k) for P in t[1]}) == 1
def zpk(t, k): return frozenset(zfree(P, k) for P in t[1])
def nbrs1(tau):
    V = list(tau); out = []
    for i in range(3):
        a, b = V[i], V[(i+1) % 3]; c = V[(i+2) % 3]
        d = (a[0] + b[0] - c[0], a[1] + b[1] - c[1])    # reflect c across edge ab
        out.append((frozenset((a, b)), frozenset((a, b, d))))
    return out
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTH[T]
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    rem = {}; creJ = {}; creA = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
        for u in s[2]: creJ.setdefault(u, j)
    for k, s in enumerate(A):
        for u in s[2]: creA.setdefault(u, k)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    bd = set(M0)
    for s in Jo[:e]: bd = (bd - s[1]) | s[2]
    Hs = Jo[e][1]; Rs = next(u for u in Hs if u[0] == 'rh'); o = min(Rs[1])
    k = ZKK[Tp]; over = {}
    for u in bd:
        if zisL(u, k): over.setdefault(zpk(u, k), []).append(u)
    assert all(len(v) == 1 for v in over.values())
    zs = [('tri', frozenset(tuple(a + b for a, b in zip(Pp, o)) for Pp in z)) for z in Zoff[Tp]]
    v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
    t3 = next(u for u in Hs if zisL(u, k) and v in u[1])
    for name, z in (('Z1', zs[0]), ('Z2', zs[1]), ('t3', t3)):
        tau = zpk(z, k)
        for (edge, nb) in nbrs1(tau):
            if nb not in over: st[('T', T, name, 'nbr outside Delta1')] += 1; continue
            u = over[nb][0]; h = tuple(int(a - b) for a, b in zip(zfix(min(u[1]), k), zfix(o, k)))
            hz = tuple(int(a - b) for a, b in zip(zfix(min(z[1]), k), zfix(o, k)))
            pv = 'H*' if u in Hs else 'creJ@%+d' % (creJ[u] - e) if u in creJ else 'creA' if u in creA else 'W' if u in W else 'X'
            hasrh = any(w[0] == 'rh' and frozenset(zfree(Pp, k) for Pp in w[1]) >= edge for w in bd)
            st[('T', T, name, 'h(self)', hz, 'nbr h', h, pv, 'rhombus over edge', hasrh)] += 1
for k in sorted(st, key=str): print(st[k], k)
