"""(Comp)(ii) cross-check on records.  Every tile removed by the re-run's middle (steps p..p+m), classified:
 R_p | common tile c (in H_p) | J_o-created | in M0 | re-run-middle-created | prefix-created.  For prefix-created: does A remove it?
 Lemma X: no tile created both by J_o and by the common prefix.  Also (T#): R_p+gamma frame cell? (re-run continued: yes by BS)
 usage: python3 comp2.py SIZE"""
import sys, pickle
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    TT = 'B' if L == 0 else 'C'; M0 = set(P['M0'])
    creJ = set().union(*[s[2] for s in Jo]); creP = set().union(*[s[2] for s in A[:p]]) if p else set()
    remA = set().union(*[s[1] for s in A])
    creMid = set()
    st[(TT, 'Lemma X: tile created by J_o and by the prefix', bool(creJ & creP))] += 1
    Rp = A[p][1]
    for i in range(p, p + len(mB)):
        for u in Bq[i][1]:
            if u in Rp: k = 'R_p' if u[0] == 'rh' else 'common tile c (in H_p)'
            elif u in creMid: k = 're-run-middle-created'
            elif u in creJ: k = 'J_o-created'
            elif u in M0: k = 'M0'
            elif u in creP: k = 'PREFIX-created, A removes it: %s' % (u in remA)
            else: k = 'OTHER'
            st[(TT, 'removed by re-run middle', u[0], k)] += 1
        creMid |= set(Bq[i][2])
for k in sorted(st, key=str): print(st[k], k)
