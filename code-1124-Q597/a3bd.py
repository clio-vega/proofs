"""(A3) instances: A-shaped wake tiles u read by A at k >= p+2 (geometric records).  For each: kind, J_o step type, and the
position of u - alpha relative to the projected triangle Delta of the lattice in which -alpha lies (k = 1 if alpha in
Lambda1 i.e. T = B, k = 0 if T = C): for a triangle, does pi(u - alpha) have an edge on the boundary of Delta (a neighbour
outside Delta)?  for a square, is its pi-edge on the boundary?   usage: python3 a3bd.py SIZE"""
import sys, os, pickle
from collections import Counter
from fractions import Fraction as Fr
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0')
def zedge(t, v): V = list(t[1]); return any(tuple(Fr(b) - Fr(a) for a, b in zip(P, Q)) == tuple(Fr(x) for x in v) for P in V for Q in V)
def zsh(t, v): return (t[0], frozenset(tuple(a + b for a, b in zip(P, v)) for P in t[1]))
def zfix(P, k): return P[:2] if k == 1 else P[2:]
def zfree(P, k): return tuple(P[2:]) if k == 1 else tuple(P[:2])
def zpk(t, k): return frozenset(zfree(P, k) for P in t[1])
def nbrs(tau):
    V = list(tau); out = []
    for i in range(3):
        a, b, c = V[i], V[(i+1) % 3], V[(i+2) % 3]; out.append(frozenset((a, b, (a[0] + b[0] - c[0], a[1] + b[1] - c[1]))))
    return out
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; k = 1 if T == 'B' else 0
    ma = (0, 0, 0, -1) if T == 'B' else (0, 1, 0, 0)
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo)); TO = typed(Jo)
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    Dl = {zpk(u, k) for u in M0 if u[0] == 'tri' and len({zfix(Q, k) for Q in u[1]}) == 1}
    for kk in range(p + 2, len(A)):
        for u in A[kk][1]:
            if u[0] == 'rh' or u not in W or not zedge(u, ma): continue
            um = zsh(u, ma)
            if um[0] == 'tri':
                tau = zpk(um, k); onb = tau in Dl and any(n not in Dl for n in nbrs(tau)); where = ('tri', 'in Delta', tau in Dl, 'on boundary', onb)
            else:
                e = zpk(um, k); where = ('sq', 'proj is edge', len(e) == 2)
                if len(e) == 2:
                    a, b = tuple(e); tris = [tau for tau in Dl if a in tau and b in tau]
                    where += ('#Delta triangles on that edge', len(tris))
            st[('T', T, 'Jo step', TO[rem[u]], 'u', u[0], 'u-alpha', where, 'in M0', um in M0)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
