"""Retreat probe (case b, (Comp)(i)): the re-run's LAST middle hexagon (= H_p + gamma under BS).  For each of its non-rh
tiles: is it in M0 & Mx?  created by J_o -- at which step relative to the split e = c0 + |mA| - 1?  created by A / the
re-run?  Also: J_o's step e is (0,b)?  Is the re-run's last middle rhombus R_{p+m+1} - alpha?   usage: python3 lastsplit.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def shift(t, v): return (t[0], frozenset(tuple(a + b for a, b in zip(P, v)) for P in t[1]))
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    rem = {}
    for j, s_ in enumerate(Jo):
        for u in s_[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); m = len(mA) - 1; e = c0 + m
    creJ = {}
    for j, s_ in enumerate(Jo):
        for u in s_[2]: creJ.setdefault(u, j)
    last = Bq[p + len(mB) - 1]
    al = disp(A[p], A[p])
    TJ = typed(Jo)
    st[('T', 'BC'[L], 'J_o step e type', TJ[e] if e < len(Jo) else None, '|mA|', len(mA), '|mB|', len(mB))] += 1
    for u in last[1]:
        if u[0] == 'rh': continue
        if u in M0 and u in Mx: kind = 'M0&Mx'
        elif u in creJ: kind = ('J_o-created at e+', creJ[u] - e)
        elif u in W: kind = 'W'
        else: kind = 'other (A/prefix/re-run created)'
        st[('T', 'BC'[L], 'last re-run middle tile', u[0], kind)] += 1
for k in sorted(st, key=str): print(st[k], k)
