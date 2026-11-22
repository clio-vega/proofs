"""case a on data: are the J_o steps read by H_p consecutive (c0, c0+1, c0+2)?  A[p+1] == rerun[p] + alpha?  rerun[p+1] == A[p] + alpha'?
suffix hypotheses at k = p+2: rhombi equal, A's later hexagons avoid D?  A reads W after p+1?   usage: python3 caseAdata.py SIZE"""
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
    if typed(Bq)[p][1] == 'b': continue
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    rem = {}
    for j, s_ in enumerate(Jo):
        for u in s_[1]: rem.setdefault(u, j)
    reads = sorted({rem[u] for u in A[p][1] if u in W})
    al = disp(A[p], A[p]); alp = disp(Bq[p], Bq[p])
    k = p + 2
    X = set(apply(M0, A[:k])); Y = set(apply(Mx, Bq[:k])); D = X - Y
    later = set().union(*[set(s_[1]) for s_ in A[k:]]) if k < len(A) else set()
    readsW = any(u in W for s_ in A[p+1:] for u in s_[1])
    readsW2 = any(u in W for s_ in A[p+2:] for u in s_[1])
    st[('consecutive', reads == list(range(reads[0], reads[0] + len(reads))), 'n', len(reads),
        'A[p+1]==rerun[p]+al', tr(Bq[p], al) == A[p+1], 'rerun[p+1]==A[p]+alp', tr(A[p], alp) == Bq[p+1],
        'R equal at p+2', (k >= len(A) and k >= len(Bq)) or (k < len(A) and k < len(Bq) and rh(A[k][1]) == rh(Bq[k][1])),
        'later avoid D', not (later & D), 'A reads W at p+1', any(u in W for u in A[p+1][1]), 'A reads W after p+1', readsW2,
        'D within W', D <= W)] += 1
for k in sorted(st, key=str): print(st[k], k)
