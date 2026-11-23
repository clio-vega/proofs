"""(Stale_a): composition of D = X_{p+2} minus Y_{p+2} in case a (L's re-run takes (1,a_T')), as stale.py does for case b.
Kinds: W (wake tile M0-Mx, not removed by A's first p+2 steps); rh; created by A at p or p+1 (non-rh) / created by the prefix;
M (tile of M0 & Mx removed by the re-run at p or p+1, not by A[:p+2]): which re-run step, and is it (tile of H_{p+1}) - alpha
or (tile of H_p) + alpha'.  Then: does any later hexagon of A (steps >= p+2) meet each kind?   usage: python3 staleA.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def shift(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
aT = ('a1', 'a0'); st = Counter(); nrec = 0
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] == 'b': continue
    nrec += 1; TT = 'B' if L == 0 else 'C'
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    k = p + 2
    X = set(apply(M0, A[:k])); Y = set(apply(Mx, Bq[:k])); D = X - Y
    al = disp(A[p], A[p]); alp = disp(Bq[p], Bq[p])
    later = set().union(*[set(s_[1]) for s_ in A[k:]]) if k < len(A) else set()
    creA = {}
    for i, s_ in enumerate(A):
        for u in s_[2]: creA.setdefault(u, i)
    st[('T', TT, '|mA|', len(mA), '|mB|', len(mB), 'len(A)-k', len(A) - k, 'len(B)-k', len(Bq) - k)] += 1
    for u in D:
        if u in W: kind = 'W'
        elif u[0] == 'rh': kind = 'rh'
        elif u in creA and creA[u] >= p: kind = ('created by A at', creA[u] - p)
        elif u in creA: kind = 'prefix-created'
        elif u in M0 and u in Mx:
            j = next((i - p for i in (p, p+1) if u in Bq[i][1]), None)
            kind = ('M: removed by rerun at p+', j, 'in H_{p+1}-alpha', u in {shift(w, tuple(-x for x in al)) for w in A[p+1][1]},
                    'in H_p+alpha\'', u in {shift(w, alp) for w in A[p][1]})
        else: kind = 'other'
        st[('T', TT, 'kind', kind, 'read later by A', u in later)] += 1
print('case-a records', nrec)
for k in sorted(st, key=str): print(st[k], k)
