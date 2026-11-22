"""Every wake tile u removed by the clean (geometric) journey A at a step k >= p+1: J_o step j(u) (class,type), k-p (capped),
local kind: B-local (u-alpha created by Jo[j(u)]), A-shaped (-alpha edge vector of u), else; u-alpha in M0; u-alpha created by
which Jo step (relative to j(u)); u-alpha in M0.Jo; u-alpha in Mx_P (WS-b).   usage: python3 wsbkind.py SIZE"""
import sys, pickle
from collections import Counter
from fractions import Fraction as Fr
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def edge_of(t, v):
    V = list(t[1]); return any(tuple(Fr(b)-Fr(a) for a, b in zip(P, Q)) == tuple(Fr(x) for x in v) for P in V for Q in V)
def shift(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx; MxP = set(apply(Mx, A[:p]))
    al = disp(A[p], A[p]); ma = tuple(-x for x in al); TO = typed(Jo)
    rem = {}; cre = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
        for u in s[2]: cre.setdefault(u, j)
    for k in range(p+1, len(A)):
        for u in A[k][1]:
            if u[0] == 'rh' or u not in W: continue
            j = rem[u]; um = shift(u, ma)
            kind = 'B-local' if um in Jo[j][2] else ('A-shaped' if edge_of(u, ma) else 'other')
            st[('T', 'B' if L == 0 else 'C', 'k-p', min(k-p, 2), 'Jo step', TO[j], 'tile', u[0], 'kind', kind,
                'u-a in M0', um in M0, 'u-a created at j-rel', None if um not in cre else cre[um] - j,
                'u-a in M0.Jo', um in Mx, 'WS-b', um in MxP)] += 1
for k in sorted(st, key=str): print(st[k], k)
