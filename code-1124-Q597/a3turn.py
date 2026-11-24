"""(Z1) = (A3) at the turn triangle.  (a) Canonical: tau3 := the triangle of H* at v that is A-shaped (-alpha an edge vector);
check tau3 - alpha == Z1, and that tau3 is the ONLY A-shaped non-rh tile of H*.  (b) Census (split records, both sizes):
the A-shaped wake tile read by A at k = p+m (J_o's (0,b) step) translated by -alpha equals the Z1 tile at the split, Z1 in M0,
Z1 never removed by J_o, Z1 not in H*.   usage: python3 a3turn.py"""
import sys, os, pickle
from collections import Counter
from fractions import Fraction as Fr
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
ZDELTA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
ZOFF = {'C': [((0, 0, 0, -1), (0, 0, 0, 0), (0, 0, 1, -1)), ((0, 0, 0, 0), (0, 0, 1, -1), (0, 0, 1, 0))],
     'B': [((0, 1, 0, 1), (0, 2, 0, 1), (1, 1, 0, 1)), ((0, 1, 0, 1), (1, 0, 0, 1), (1, 1, 0, 1))]}
def zedge_of(t, v): V = list(t[1]); return any(tuple(Fr(b) - Fr(a) for a, b in zip(P, Q)) == tuple(Fr(x) for x in v) for P in V for Q in V)
def zsh(t, v): return (t[0], frozenset(tuple(Fr(a) + Fr(b) for a, b in zip(P, v)) for P in t[1]))
for T in 'BC':
    Tp = OTHER[T]; ma = tuple(-x for x in ZDELTA[T])
    Rs = CF[Tp][(k0, 'b')][0]; HJ, RJ, CJ = step(Tp, Rs, 'b'); o = min(Rs[1])
    z1 = ('tri', frozenset(tuple(Fr(a) + b for a, b in zip(P, o)) for P in ZOFF[Tp][0]))
    ash = [u for u in HJ if u[0] != 'rh' and zedge_of(u, ma)]
    print('(a) T', T, '| A-shaped non-rh tiles of H*:', [u[0] for u in ash], '| tau3 - alpha == Z1:', [zsh(u, ma) == z1 for u in ash])
sys.argv = ['x', '5', '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
aT = ('a1', 'a0'); st = Counter()
for SIZE in (5, 6):
    PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
    for P in PAIRS:
        if P['kind'] != 'BC': continue
        T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
        L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
        t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
        if typed(Bq)[p][1] != 'b': continue
        T = 'BC'[L]; Tp = OTHER[T]; ma = ZDELTA[T]; ma = tuple(-x for x in ma)
        M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
        rem = {}
        for j, s in enumerate(Jo):
            for u in s[1]: rem.setdefault(u, j)
        c0 = min(rem[u] for u in A[p][1] if u in W); m = len(mA) - 1; e = c0 + m
        Rs = next(u for u in Jo[e][1] if u[0] == 'rh'); o = min(Rs[1])
        z1 = ('tri', frozenset(tuple(a + b for a, b in zip(P_, o)) for P_ in ZOFF[Tp][0]))
        ash = [u for u in A[p + m][1] if u[0] != 'rh' and u in W and zedge_of(u, ma)]
        ok = [tuple(tuple(int(x) for x in Pp) for Pp in sorted(zsh(u, ma)[1])) == tuple(sorted(z1[1])) for u in ash]
        st[('size', SIZE, 'T', T, 'm>=2', m >= 2, '#A-shaped read at p+m', len(ash), 'u-alpha == Z1', tuple(ok), 'Z1 in M0', z1 in M0,
            'J_o removes Z1', z1 in rem, 'Z1 in H*', z1 in Jo[e][1])] += 1
for k in sorted(st, key=str): print(st[k], k)
