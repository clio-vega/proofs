"""PROVE 11-17 c2 step 1.  For each pair (x->B, y->C) [and B/B control x<y] on the Q594 sampler (census, states cap 60),
compute pure-cell data (q,p) of M, M.x, M.y, M.x.y, M.y.x and the NET change of each move as a dict
(i,t) -> (old vertex, new vertex), i=0: q (t in Tri(P0), value in Vert(P1)); i=1: p.
Pickles records for analysis.  Exact integers only."""
import sys, copy, pickle
sys.path.insert(0, '/home/clio/projects/proofs/code-1116-Q585-MC')
from qp import *          # brings rec.py -> common.py (census, corners, mig1, tiles_of, states, mig_rec)
def pq(b):
    q, p = qp(b); return {**{(0, t): v for t, v in q.items()}, **{(1, s): v for s, v in p.items()}}
def delta(A, B):
    assert set(A) == set(B)
    return {k: (A[k], B[k]) for k in A if A[k] != B[k]}
def mv(b, x, T, DV):
    bb = copy.deepcopy(b); mig1(bb, 0, x, DV[T]); return bb
recs = []
for key, b0 in census():
    for b, DV in states(b0, cap=60):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind, (tx, ty) in (('BC', 'BC'), ('BB', 'BB')):
                    if kind == 'BB' and not x < y: continue
                    b1 = mv(b, x, tx, DV); b2 = mv(b, y, ty, DV)
                    b12 = mv(b1, y, ty, DV); b21 = mv(b2, x, tx, DV)
                    P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                    com = tiles_of(b12) == tiles_of(b21)
                    recs.append(dict(kind=kind, key=key, x=x, y=y, com=com, pqsame=(P12 == P21),
                                     X=delta(P, P1), Y=delta(P, P2), Yx=delta(P1, P12), Xy=delta(P2, P21),
                                     P=P))
pickle.dump(recs, open('recs.pkl', 'wb'))
from collections import Counter
print(Counter((r['kind'], r['com'], r['pqsame']) for r in recs))
