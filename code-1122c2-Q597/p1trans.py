"""Step p+1 of L (the (1,b) step) translated by -alpha equals the re-run's step p (geometric records). usage: python3 p1trans.py SIZE"""
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
    al = disp(A[p], A[p]); ma = tuple(-x for x in al)
    st[('A[p+1] type', typed(A)[p+1], 'tr(A[p+1], -alpha) == rerun[p]', tr(A[p+1], ma) == Bq[p])] += 1
    # lesion: translate by +alpha instead
    st[('LESION tr(A[p+1], +alpha) == rerun[p]', tr(A[p+1], al) == Bq[p])] += 1
for k in sorted(st, key=str): print(st[k], k)
