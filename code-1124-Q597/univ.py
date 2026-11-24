"""Tiles near R* on J_o's board right after its (0,b) step (board at e+1), in canonical position (tile translated by
-min R*).  For each target T' of J_o: positions present in ALL split records vs frequency at generic (0,b) steps followed
by the straight step.  Report universal-at-split tiles and whether they conflict with the straight configuration S.
Also provenance at the split (M0&W / M0&Mx / J_o-created step rel e).   usage: python3 univ.py SIZE"""
import sys, os, pickle
from collections import Counter, defaultdict
SIZE = int(sys.argv[1])
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import conflict, phi, poly
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
STRT = {'B': 'a1', 'C': 'a0'}; aT = ('a1', 'a0')
def nz(t): return (t[0], frozenset(tuple(Fr(x) for x in P) for P in t[1]))
def sh(t, v): return (t[0], frozenset(tuple(Fr(a) - b for a, b in zip(P, v)) for P in t[1]))
cnt = {k: Counter() for k in ('split', 'gstr')}; N = Counter(); prov = defaultdict(Counter)
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
            else: continue
            Rs = next(u for u in s[1] if u[0] == 'rh'); o = min(nz(Rs)[1]); c = phi(o)
            N[(lab, Tp)] += 1
            for u in bd:
                if min(abs(z - c) for z in poly(nz(u))) > 2.6: continue
                key = sh(nz(u), o); cnt[lab][(Tp, key)] += 1
                if lab == 'split':
                    prov[(Tp, key)]['M0&W' if u in W else 'M0&Mx' if u in M0 else 'creJ@e%+d' % (creJ.get(u, -99) - j)] += 1
print('records', dict(N))
for Tp in 'CB':
    Rs = CF[Tp][(k0, 'b')][0]; HJ, RJ, CJ = step(Tp, Rs, 'b'); S, _, _ = step(Tp, RJ, STRT[Tp])
    o = min(Rs[1])
    U = [k for (tp, k), n in cnt['split'].items() if tp == Tp and n == N[('split', Tp)]]
    print('==== J_o toward', Tp, ': universal-at-split tiles', len(U))
    for k in sorted(U, key=lambda k: (k[0], sorted(k[1]))):
        kk = (k[0], frozenset(tuple(a + b for a, b in zip(P, o)) for P in k[1]))
        g = cnt['gstr'][(Tp, k)]
        print('  ', k[0], sorted(tuple(int(x) for x in P) for P in k[1]), '| generic-straight freq %d/%d' % (g, N[('gstr', Tp)]),
              '| conflicts S:', bool(conflict(S, {kk})), '| in CJ:', kk in CJ, '| prov', dict(prov[(Tp, k)]))
