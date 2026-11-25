"""type word of A's middle a.C (steps p..p+m) and of the next step, on geometric records.  usage: python3 midword.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    TA = typed(A); m = len(mA) - 1
    w = ws(TA[p:p+m+1]); nxt = ws(TA[p+m+1:p+m+2]) or 'END'
    import re
    shape = re.sub(r'(0a[01] )+', '(0a)* ', w + ' ')
    st[('B' if L == 0 else 'C', 'middle word shape', shape.strip(), 'next', nxt)] += 1
    st[('B' if L == 0 else 'C', 'middle word', w)] += 1
for k in sorted(st, key=str): print(st[k], k)
