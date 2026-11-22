"""WS_last probe: for geometric records, where are the tiles of H_p translated by gamma (gamma = disp(C[0], C[-1]))?
For each non-rh tile t of H_p: kind (W/common), t+gamma: in M0? created by J_o at which step relative to e (J_o's turn)?
created by A? on Y_{p+m} (re-run board before its last middle step)?   usage: python3 wslast.py SIZE"""
import sys, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def shift(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx; TO = typed(Jo)
    gam = disp(mA[1], mA[-1]); k = p + len(mA) - 1
    Yk = set(apply(Mx, Bq[:k]))
    creJ = {}; remJ = {}
    for j, s_ in enumerate(Jo):
        for u in s_[2]: creJ.setdefault(u, j)
        for u in s_[1]: remJ.setdefault(u, j)
    # J_o's turn e = J_o step removing the last shadow hexagon = A[k-?]; take e = remJ of a W tile of A[p+len(mA)-1]
    eW = [remJ[u] for u in A[k][1] if u in W]
    e = max(eW) if eW else None
    for u in A[p][1]:
        if u[0] == 'rh': continue
        ug = shift(u, gam)
        where = 'M0' if ug in M0 else ('Jo-created at e+%s' % (creJ[ug] - e) if ug in creJ and e is not None else ('Jo-created' if ug in creJ else 'neither'))
        st[('T', 'B' if L == 0 else 'C', 'tile', u[0], 'in W', u in W, 't+gamma', where, 'removed by Jo later', ug in remJ and (ug not in creJ or remJ[ug] > creJ[ug]), 'on Y_{p+m}', ug in Yk)] += 1
for k in sorted(st, key=str): print(st[k], k)
