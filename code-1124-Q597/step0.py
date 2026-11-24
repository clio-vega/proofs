"""Step 0 (PROVE 11-24): is (E1) the same statement as (Tn)?  For every case-b record (L's re-run takes b):
 e_E1 = c0 + |mA| - 1 (lastsplit.py: J_o's image of the split rhombus), (E1) := J_o has step e_E1+1 and it is not (1,a_T').
 e_Tn = the J_o steps that remove a wake SQUARE read by A after p and are (0,b) (tnext.py), (Tn) := J_o has step e+1,
        not J_o's straight step (1,a_T').  Prints the 2x2 table and the index agreement.   usage: python3 step0.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0'); STR = {0: (1, 'a0'), 1: (1, 'a1')}   # straight step of J_o: L=0 (T=B) -> J_o toward C -> (1,a0)
st = Counter(); n = 0
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    n += 1
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo)); TO = typed(Jo)
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); eE = c0 + len(mA) - 1
    es = sorted({rem[u] for k in range(p+1, len(A)) for u in A[k][1] if u[0] == 'sq' and u in W and TO[rem[u]] == (0, 'b')})
    holds = lambda e: e + 1 < len(Jo) and TO[e+1] != STR[L]
    E1 = holds(eE)
    st[('#es', len(es))] += 1
    st[('index agree', es == [eE])] += 1
    st[('J_o step eE type', TO[eE] if eE < len(Jo) else None)] += 1
    st[('(0,b) steps of J_o total', sum(1 for x in TO if x == (0, 'b')))] += 1
    for e in es: st[('2x2', 'Tn ' + str(holds(e)), 'E1 ' + str(E1))] += 1
    if es != [eE]: print('DIFFER', P.get('id'), L, eE, es, TO)
print('case-b records', n)
for k in sorted(st, key=str): print(st[k], k)
