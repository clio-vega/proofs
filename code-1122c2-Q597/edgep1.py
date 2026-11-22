"""Is -alpha an edge vector of L's rhombus R_k for k = p+1 (class 1) and for k = p+2..p+|mA|-? (L-records, re-run b).
usage: python3 edgep1.py SIZE"""
import sys, pickle
from collections import Counter
from fractions import Fraction as Fr
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def edgevecs(R):
    V = list(R); out = set()
    for P in V:
        for Q in V:
            d = tuple(Fr(b)-Fr(a) for a, b in zip(P, Q))
            if any(d) and (any(d[:2]) != any(d[2:])): out.add(d)
    return out
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = []
    for i, (t, A, Bq, Jo) in enumerate(P['roles']): T_.append(typed(A)[middles(A, Bq)[2]])
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    al = disp(A[p], A[p]); ma = tuple(-x for x in al)
    TA = typed(A)
    for k in range(p+1, p+len(mA)):
        st[('T', 'B' if L == 0 else 'C', 'k-p', min(k-p, 2), 'class', TA[k][0], '-alpha edge of R_k', ma in edgevecs(rh(A[k][1])))] += 1
for k in sorted(st, key=str): print(st[k], k)
