"""(Stale): composition of D minus W at the rejoin k = p+|mA| (geometric records): each tile classified as
rh / created by A's middle (non-rh) / prefix-created / M0&Mx removed by the re-run middle: at which re-run middle index
(0..|mB|-1) and is it  (shadow tile of A) - alpha ?   usage: python3 stale.py SIZE"""
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
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    k = p + len(mA)
    X = set(apply(M0, A[:k])); Y = set(apply(Mx, Bq[:k])); D = X - Y
    al = disp(A[p], A[p]); ma = tuple(-x for x in al)
    shadowtiles = set().union(*[set(A[i][1]) for i in range(p+2, k)]) if k > p+2 else set()
    creA = {}
    for i, s_ in enumerate(A):
        for u in s_[2]: creA.setdefault(u, i)
    remB = {}
    for i, s_ in enumerate(Bq):
        for u in s_[1]: remB.setdefault(u, i)
    for u in D - W:
        if u[0] == 'rh': kind = 'rh'
        elif u in creA and creA[u] >= p: kind = 'created by A middle'
        elif u in creA: kind = 'prefix-created'
        elif u in M0 and u in Mx and u in remB:
            j = remB[u] - p
            kind = ('M: removed by rerun middle idx', 'last' if j == len(mB) - 1 else ('first' if j == 0 else 'mid'),
                    'is (A shadow tile - alpha)', any(shift(w, ma) == u for w in shadowtiles),
                    'is (A[p] tile + gamma)', u in {shift(w, disp(mA[1], mA[-1])) for w in A[p][1]} if len(mA) > 1 else None)
        else: kind = 'other'
        st[('T', 'B' if L == 0 else 'C', kind)] += 1
for k in sorted(st, key=str): print(st[k], k)
