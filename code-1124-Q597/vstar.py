"""Star of the projected vertex vbar = proj(v) in Delta (Delta1 for J_o->C, Delta0 for J_o->B) on J_o's board before step e,
split records.  For each unit triangle of Delta around vbar (ccw from tau3 = H*'s triangle): height - h(tau3), provenance;
for each edge at vbar: the stack over it (#squares, #rhombi, and which rhombi: R*, nest).  Also: is vbar on the boundary
of Delta (fewer than 6 triangles)?   usage: python3 vstar.py SIZE"""
import sys, os, pickle, math, cmath
from collections import Counter
SIZE = int(sys.argv[1])
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0'); OTH = {'B': 'C', 'C': 'B'}; ZKK = {'C': 1, 'B': 0}
def zfix(P, k): return P[:2] if k == 1 else P[2:]
def zfree(P, k): return P[2:] if k == 1 else P[:2]
def zisL(t, k): return t[0] == 'tri' and len({zfix(P, k) for P in t[1]}) == 1
def zpk(t, k): return frozenset(zfree(P, k) for P in t[1])
W6 = cmath.exp(1j * math.pi / 3)
def ang(a, b): z = complex(b[0] - a[0], 0) + (b[1] - a[1]) * W6; return math.degrees(cmath.phase(z)) % 360
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
    vb = zfree(v, k); h0 = zfix(v, k)
    over = {zpk(u, k): u for u in bd if zisL(u, k)}
    t3 = next(u for u in Hs if zisL(u, k) and v in u[1]); a3 = None
    tris = []
    for tau, u in over.items():
        if vb in tau:
            c = [x for x in tau if x != vb]; mid = ((c[0][0] + c[1][0]) / 2, (c[0][1] + c[1][1]) / 2)
            hh = tuple(int(a - b) for a, b in zip(zfix(min(u[1]), k), h0))
            pv = 'H*' if u in Hs else 'creJ@%+d' % (creJ[u] - e) if u in creJ else 'W' if u in W else 'X'
            tris.append((ang(vb, mid), hh, pv))
    a3 = next(a for a, hh, pv in tris if pv == 'H*')
    tris = sorted(((round((a - a3) % 360), hh, pv) for a, hh, pv in tris))
    # stacks over edges at vbar
    stacks = []
    for d in [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]:
        nb = (vb[0] + d[0], vb[1] + d[1]); edge = frozenset((vb, nb))
        sq = sum(1 for u in bd if u[0] == 'sq' and zpk(u, k) == edge)
        rhs = ['R*' if u == Rs else 'rh' for u in bd if u[0] == 'rh' and zpk(u, k) == edge]
        if sq or rhs: stacks.append((round((ang(vb, nb) - a3) % 360), sq, tuple(rhs)))
    st[('T', T, '#tris at vbar', len(tris), 'tris(angle from tau3, h, prov)', tuple(tris), 'stacks(angle,#sq,rh)', tuple(sorted(stacks)))] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
