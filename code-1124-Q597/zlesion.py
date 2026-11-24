"""LESION for the reduction (E1) <= (Z): drop 'this is a split' and test on EVERY (0,b) step of every J_o (both roles):
is the triangle pair Z = {Z1, Z2} (canonical offsets from min R*, from univ5.out) on J_o's board at e+1, and what is the
next step?  Expect: Z present => next != straight (geometric, Z conflicts S); generic straight steps have Z absent.
usage: python3 zlesion.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1])
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
STRT = {'B': (1, 'a1'), 'C': (1, 'a0')}; aT = ('a1', 'a0')
Z = {'C': [((0, 0, 0, -1), (0, 0, 0, 0), (0, 0, 1, -1)), ((0, 0, 0, 0), (0, 0, 1, -1), (0, 0, 1, 0))],
     'B': [((0, 1, 0, 1), (0, 2, 0, 1), (1, 1, 0, 1)), ((0, 1, 0, 1), (1, 0, 0, 1), (1, 1, 0, 1))]}
def nz(t): return (t[0], frozenset(tuple(Fr(x) for x in P) for P in t[1]))
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    for role in (0, 1):
        t, A, Bq, Jo = P['roles'][role]; T = 'BC'[role]; Tp = OTHER[T]; TO = typed(Jo)
        mA, mB, p, q = middles(A, Bq); M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
        split = None
        if role == L and typed(Bq)[p][1] == 'b':
            rem = {}
            for j, s in enumerate(Jo):
                for u in s[1]: rem.setdefault(u, j)
            split = min(rem[u] for u in A[p][1] if u in W) + len(mA) - 1
        bd = set(M0)
        for j, s in enumerate(Jo):
            bd = (bd - s[1]) | s[2]
            if TO[j] != (0, 'b'): continue
            o = min(nz(next(u for u in s[1] if u[0] == 'rh'))[1])
            nb = {nz(u) for u in bd}
            zs = [('tri', frozenset(tuple(a + b for a, b in zip(P_, o)) for P_ in z)) for z in Z[Tp]]
            pres = tuple(z in nb for z in zs)
            nxt = 'END' if j + 1 == len(Jo) else ('STRAIGHT' if TO[j+1] == STRT[Tp] else 'other')
            st[('split' if j == split else 'generic', 'J_o->' + Tp, 'Z1,Z2 present', pres, 'next', nxt)] += 1
for k in sorted(st, key=str): print(st[k], k)
