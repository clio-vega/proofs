"""PROVE 11-25: WS_last probe 2. For each geometric record (role L, re-run step p of type b), each non-rh tile tau of H_p,
t = tau + gamma (the re-run's last-middle hexagon tile) and t + alpha (the same tile at A's side).
Report: is t removed by a C-alpha step of the re-run (want 0), and the status of t+alpha w.r.t. A: in M0? removed by A at which
step (relative to p)? created by A at which step? and w.r.t. J_o.   usage: python3 wsl1.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
def shift(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
aT = ('a1', 'a0'); st = Counter(); nrec = 0
def idx(J):
    cre = {}; rem = {}
    for j, s_ in enumerate(J):
        for u in s_[2]: cre.setdefault(u, j)
        for u in s_[1]: rem.setdefault(u, j)
    return cre, rem
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    nrec += 1
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    al = disp(mA[0], mA[0]); gam = disp(mA[1], mA[-1]); m = len(mA) - 1
    creA, remA = idx(A); creB, remB = idx(Bq); creJ, remJ = idx(Jo)
    TT = 'B' if L == 0 else 'C'
    # sanity: middle of re-run is C - alpha then a + gamma
    ok = all(Bq[p+i-1][1] == frozenset(shift(u, tuple(-x for x in al)) for u in A[p+i][1]) for i in range(1, m+1))
    st[(TT, 'pre: re-run steps p..p+m-1 = A steps p+1..p+m translated by -alpha', ok)] += 1
    st[(TT, 'pre: m', m)] += 1
    for tau in A[p][1]:
        if tau[0] == 'rh': continue
        tg = shift(tau, gam); ta = shift(tg, al)
        hit = [i for i in range(1, m+1) if tg in Bq[p+i-1][1]]
        st[(TT, 'WS_last: t removed by a C-alpha step', bool(hit))] += 1
        s = []
        s.append('t+a in M0' if ta in M0 else 't+a notin M0')
        s.append('remA@%s' % (remA[ta]-p if ta in remA else None))
        s.append('creA@%s' % (creA[ta]-p if ta in creA else None))
        s.append('remJ@%s' % (remJ[ta] if ta in remJ else None))
        s.append('creJ' if ta in creJ else '-')
        st[(TT, tau[0], 'tau in W', tau in W, 'm', 'any', *s)] += 1
print('records', nrec)
for k in sorted(st, key=str): print(st[k], k)
