"""Check the first-moment proof of (N) on data: for L with re-run b, M(Delta) = sum_k pos(k) (x) Delta(k) on side s
(pos = lexmin vertex of the key triangle), entry [1][1] predicted = -(#1b + 2#0a1 + #0b) (T=B) / -(2#0a0 + #0b + #1b) (T=C).
Also: supp Delta within side(L) in BOTH cases.   usage: python3 momcheck.py SIZE"""
import sys, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def sub(D, E):
    out = dict(D)
    for k, v in E.items():
        w = tuple(a-b for a, b in zip(out.get(k, (0, 0)), v)); out[k] = w
        if w == (0, 0): del out[k]
    return out
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    Delta = sub(P['Xy'], P['X'])
    T_ = []
    for i, (t, A, Bq, Jo) in enumerate(P['roles']): T_.append(typed(A)[middles(A, Bq)[2]])
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    s = 1 if L == 0 else 0
    case = typed(Bq)[p][1]
    st[('supp Delta within side(L)', case, all(k[0] == s for k in Delta))] += 1
    if case != 'b': continue
    M = [[0, 0], [0, 0]]
    for (side, key), val in Delta.items():
        pos = min(key)
        for i in range(2):
            for j in range(2): M[i][j] += pos[i] * val[j]
    TC = Counter(typed(A)[p+1:p+len(mA)])
    pred = -(TC[(1, 'b')] + 2*TC[(0, 'a1')] + TC[(0, 'b')]) if L == 0 else -(2*TC[(0, 'a0')] + TC[(0, 'b')] + TC[(1, 'b')])
    st[('case b', 'T', 'B' if L == 0 else 'C', 'M[1][1] == predicted', M[1][1] == pred, 'M != 0', any(any(r) for r in M), 'pred<=-1', pred <= -1)] += 1
for k in sorted(st, key=str): print(st[k], k)
