"""For split records: m distribution; for the two universal triangles Z1 (M0&Mx) and Z2 (W) [the tiles of J_o's board at
e+1 that conflict with the straight configuration S]: which steps of A / J_o have a hexagon or created set sharing an EDGE
with them (relative to p+m resp. e), and when J_o removes Z2.  Also: does S conflict with A's tiles created at ANY step?
usage: python3 zadj.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1])
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import conflict
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
STRT = {'B': 'a1', 'C': 'a0'}; aT = ('a1', 'a0')
def nz(t): return (t[0], frozenset(tuple(Fr(x) for x in P) for P in t[1]))
def edges(t):
    V = list(t[1]); return {frozenset((a, b)) for a in V for b in V if a < b and sum(abs(float(x) - float(y)) for x, y in zip(a, b)) <= 2.01 and
                           ((a[0] == b[0] and a[1] == b[1]) or (a[2] == b[2] and a[3] == b[3]))}
def adj(u, S): return any(edges(u) & edges(s) for s in S if s != u)
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTHER[T]; TO = typed(Jo)
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); m = len(mA) - 1; e = c0 + m
    st[('T', T, 'm', m)] += 1
    bd = set(M0)
    for s in Jo[:e+1]: bd = (bd - s[1]) | s[2]
    RJ = next(u for u in Jo[e][2] if u[0] == 'rh')
    S, _, _ = step(Tp, nz(RJ), STRT[Tp])
    nb = {nz(u): u for u in bd}
    Z = {nb[b] for (_, b) in conflict(S, set(nb)) if b[0] == 'tri' and nb[b] in M0}
    for u in Z:
        kind = 'Z1(Mx)' if u in Mx else 'Z2(W)@e%+d' % (rem[u] - e)
        nu = nz(u)
        aA = tuple(sorted(k - p - m for k, s in enumerate(A) if adj(nu, {nz(x) for x in s[1] | s[2]})))
        aJ = tuple(sorted(j - e for j, s in enumerate(Jo) if adj(nu, {nz(x) for x in s[1] | s[2]})))
        st[('T', T, kind, 'A steps adj (rel p+m)', aA[:6], 'Jo steps adj (rel e)', aJ[:6])] += 1
    creA = set()
    for s in A: creA |= {nz(x) for x in s[2]}
    st[('T', T, 'S conflicts some A-created tile', bool(conflict(S, creA)), 'A ends at R#', p + len(mA) == len(A))] += 1
for k in sorted(st, key=str): print(st[k], k)
