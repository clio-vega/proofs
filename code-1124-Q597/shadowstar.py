"""Along the shadow j in [c0+2, e] of split records: J_o's step-j type, and the vertex star at v_j (v = min R^o_j for J_o->C,
min+(0,1,0,1) for J_o->B) on J_o's board BEFORE step j, as (sector start, width, kind, provenance) with provenance
X = M0&Mx, W(+d) = removed by J_o at j+d, j(-d) = J_o-created at j-d.   usage: python3 shadowstar.py SIZE"""
import sys, os, pickle, math, cmath
from collections import Counter
SIZE = int(sys.argv[1])
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import phi
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0')
def nz(t): return (t[0], frozenset(tuple(Fr(x) for x in P) for P in t[1]))
def sector(t, v):
    zv = phi(v); pts = [phi(P) for P in t[1]]
    if not any(abs(z - zv) < 1e-9 for z in pts): return None
    a = sorted(math.degrees(cmath.phase(z - zv)) % 360 for z in pts if abs(abs(z - zv) - 1) < 1e-6)
    cdir = math.degrees(cmath.phase(sum(pts) / len(pts) - zv)) % 360
    for s, e_ in ((a[0], a[1]), (a[1], a[0])):
        w = (e_ - s) % 360
        if (cdir - s) % 360 < w: return (round(s) % 360, round(w))
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTHER[T]; TO = typed(Jo)
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    rem = {}; creJ = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
        for u in s[2]: creJ.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); m = len(mA) - 1; e = c0 + m
    bd = set(M0)
    for j, s in enumerate(Jo):
        if c0 + 2 <= j <= e:
            Rs = next(u for u in s[1] if u[0] == 'rh'); o = min(nz(Rs)[1])
            v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
            star = []
            for u in bd:
                sc = sector(nz(u), v)
                if sc is None: continue
                pv = 'X' if u in Mx and u in M0 else 'W%+d' % (rem[u] - j) if u in W else 'j%+d' % (creJ.get(u, -99) - j)
                star.append((sc[0], sc[1], u[0], pv))
            star.sort()
            st[('T', T, 'pos', 'e' if j == e else 'shadow<e', TO[j], tuple(star))] += 1
        bd = (bd - s[1]) | s[2]
for k in sorted(st, key=str): print(st[k], k)
