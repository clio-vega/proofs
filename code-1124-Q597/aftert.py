"""After the turn: does A take steps beyond p+m, and does it read wake tiles there (any shape)?  usage: python3 aftert.py SIZE"""
import sys
from collections import Counter
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/bandrule.py').read().split('st = Counter(); shapes')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1; m = e - c0
    after = len(A) - 1 - (p + m)
    wr = sum(1 for kk in range(p + m + 1, len(A)) for u in A[kk][1] if u in W and u[0] != 'rh')
    st[(T, 'A steps after p+m', min(after, 3), 'wake non-rh reads after', min(wr, 3))] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
