"""Which tile of J_o's board at step e+1 blocks J_o's straight configuration at R^o_{e+1}?  Case-b records.
Straight configuration placed by words.configs (hsderive preamble).  For each blocker: offset of its lexmin from min R*,
kind, provenance (M0 / J_o-created at step j-e / in A's board at p+m / removed by A at step k-p / nest rhombus).
CONTROL: the same at the generic (0,b) steps of all J_o's that ARE followed by the straight step: 0 blockers expected.
usage: python3 blocker.py SIZE"""
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
def sub(P, Q): return tuple(int(a - b) for a, b in zip(P, Q))
st = Counter(); ctl = Counter(); EX = []
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    T = 'BC'[L]; Tp = OTHER[T]
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx; TO = typed(Jo)
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    creJ = {}
    for j, s in enumerate(Jo):
        for u in s[2]: creJ.setdefault(u, j)
    remA = {}
    for k, s in enumerate(A):
        for u in s[1]: remA.setdefault(u, k)
    caseb = typed(Bq)[p][1] == 'b'
    c0 = min(rem[u] for u in A[p][1] if u in W) if caseb else None
    bd = set(M0)
    for j, s in enumerate(Jo):
        bd = (bd - s[1]) | s[2]
        if TO[j] != (0, 'b') or j + 1 >= len(Jo): continue
        RJ = next(u for u in s[2] if u[0] == 'rh')
        Rstar = next(u for u in s[1] if u[0] == 'rh')
        S, _, _ = step(Tp, nz(RJ), STRT[Tp])
        nb = {nz(u): u for u in bd}
        blk = conflict(S, set(nb))
        issplit = caseb and j == c0 + len(mA) - 1
        if issplit:
            st[('T', T, 'split: #blockers', len(blk), 'A ends at R#', p + len(mA) == len(A))] += 1
            for (sx, b) in blk:
                ob = nb[b]; off = sub(min(b[1]), min(nz(Rstar)[1]))
                prov = []
                if ob in M0: prov.append('M0' + ('&W' if ob in W else '&Mx'))
                if ob in creJ: prov.append('creJ@e%+d' % (creJ[ob] - j))
                if ob in remA: prov.append('remA@p%+d' % (remA[ob] - p))
                st[('T', T, 'blocker', b[0], off, 'vs straight tile', sx[0], tuple(prov))] += 1
        else:
            ctl[('T', Tp, 'generic (0,b), next=straight' if TO[j+1] == (1, STRT[Tp]) else 'generic, next=other', 'blocked', bool(blk))] += 1
for k in sorted(st, key=str): print(st[k], k)
print('--- CONTROL (non-split (0,b) steps of J_o)')
for k in sorted(ctl, key=str): print(ctl[k], k)
