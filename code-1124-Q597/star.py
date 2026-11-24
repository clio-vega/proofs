"""Full vertex star at the hexagon vertex v (v = o for J_o->C, o+(0,1,0,1) for J_o->B) on J_o's board at e+1, split
records vs generic (0,b)->straight; each tile labelled (kind, sector start, width, provenance).  Also A's involvement:
for split records, whether A's board at p+m+1 has the same star outside the hexagon region.   usage: python3 star.py SIZE"""
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
STRT = {'B': 'a1', 'C': 'a0'}; aT = ('a1', 'a0')
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
    for role in (0, 1):
        t, A, Bq, Jo = P['roles'][role]; T = 'BC'[role]; Tp = OTHER[T]; TO = typed(Jo)
        mA, mB, p, q = middles(A, Bq)
        M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
        split = None
        if role == L and typed(Bq)[p][1] == 'b':
            rem = {}
            for j, s in enumerate(Jo):
                for u in s[1]: rem.setdefault(u, j)
            split = min(rem[u] for u in A[p][1] if u in W) + len(mA) - 1
        creJ = {}
        for j, s in enumerate(Jo):
            for u in s[2]: creJ.setdefault(u, j)
        bd = set(M0)
        for j, s in enumerate(Jo):
            bd = (bd - s[1]) | s[2]
            if TO[j] != (0, 'b') or j + 1 >= len(Jo): continue
            if j == split: lab = 'split'
            elif TO[j+1] == (1, STRT[Tp]): lab = 'gstr'
            else: lab = 'gother'
            Rs = next(u for u in s[1] if u[0] == 'rh'); o = min(nz(Rs)[1])
            v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
            star = []
            for u in bd:
                sc = sector(nz(u), v)
                if sc is None: continue
                pv = 'W' if u in W else 'X' if u in M0 else 'j%+d' % (creJ.get(u, -99) - j)
                star.append((sc[0], sc[1], u[0], pv))
            star.sort()
            st[(lab, Tp, tuple(star))] += 1
for k in sorted(st, key=str): print(st[k], k)
