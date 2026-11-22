"""WS-a on the shadow range: for geometric records, every non-rh tile of H_k, p+1 <= k <= p+|mA|-1, is in W or created by A at a step
in [p, k).  Also: tiles created by A's prefix ever read at k >= p+1?  usage: python3 wsa.py SIZE"""
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
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
    creA = {}
    for i, s_ in enumerate(A):
        for u in s_[2]: creA.setdefault(u, i)
    for k in range(p+1, len(A)):
        rng = 'shadow' if k <= p + len(mA) - 1 else 'after middle'
        for u in A[k][1]:
            if u[0] == 'rh': continue
            c = creA.get(u)
            kind = 'W' if u in W else ('middle-made' if c is not None and p <= c < k else ('prefix-made' if c is not None and c < p else 'common M0 (not W)'))
            st[(rng, kind)] += 1
for k in sorted(st, key=str): print(st[k], k)
