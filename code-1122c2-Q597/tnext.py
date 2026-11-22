"""J_o's step after the turn b-step e (the class-0 b-step whose square A reads), geometric records. usage: python3 tnext.py SIZE"""
import sys, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo)); TO = typed(Jo)
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    es = sorted({rem[u] for k in range(p+1, len(A)) for u in A[k][1] if u[0] == 'sq' and u in W and TO[rem[u]] == (0, 'b')})
    for e in es:
        st[('T', 'B' if L == 0 else 'C', 'Jo after turn', TO[e+1] if e+1 < len(Jo) else 'END')] += 1
for k in sorted(st, key=str): print(st[k], k)
