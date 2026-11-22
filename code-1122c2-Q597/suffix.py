"""Suffix Lemma dinv hypothesis (MR) for L-records with re-run b: at k = first common-suffix step (= p+|mA|),
R_k = R'_k; D = X_k minus Y_k, D' = Y_k minus X_k; does any later hexagon of A meet D (or of A' meet D')?
Also: D within W (wake tiles not removed by A) / within tiles created by A's middle; and the same at k-1 (lesion: one step early).
usage: python3 suffix.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    M0 = set(P['M0']); T = []
    for i, (t, A, Bq, Jo) in enumerate(P['roles']): T.append(typed(A)[middles(A, Bq)[2]])
    L = [i for i in (0, 1) if T[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]
    mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    Mx = set(apply(M0, Jo)); W = M0 - Mx
    for lab, k in (('k=p+|mA|', p + len(mA)), ('LESION k-1', p + len(mA) - 1)):
        X = set(apply(M0, A[:k])); Y = set(apply(Mx, Bq[:k]))
        D, D2 = X - Y, Y - X
        sameR = k < len(A) and k < len(Bq) and rh(A[k][1]) == rh(Bq[k][1]) if k < min(len(A), len(Bq)) else ('end', len(A) == k, len(Bq) == k)
        laterA = set().union(*[set(s[1]) for s in A[k:]]) if k < len(A) else set()
        laterB = set().union(*[set(s[1]) for s in Bq[k:]]) if k < len(Bq) else set()
        createdA = set().union(*[set(s[2]) for s in A[p:k]]) if k > p else set()
        st[(lab, 'R_k=R\'_k', sameR, 'A later meets D', bool(laterA & D), "A' later meets D'", bool(laterB & D2),
            'D in W', D <= W, 'D minus W in created-by-A-middle', (D - W) <= createdA)] += 1
for k in sorted(st, key=str): print(st[k], k)
