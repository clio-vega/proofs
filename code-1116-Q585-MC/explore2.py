"""For each diamond (x->B, y->C) [and control (x->B, y->B)] compare the displacement of each move in the two orders.
D(move, state) = (dq, dp) dicts.  Tests: moved SETS equal across orders?  moved VECTORS equal?  and additivity (=MC)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1115c2-Q585-G1'); sys.path.insert(0,'.')
from rec import *
from explore1 import qp
from collections import Counter
def disp(b, x, dvec):
    q0, p0 = qp(b); Y, ij, h = mig_rec(b, 0, x, dvec); q1, p1 = qp(b)
    dq = {t: tuple(a-c for a, c in zip(q1[t], q0[t])) for t in q0 if q1[t] != q0[t]}
    dp = {s: tuple(a-c for a, c in zip(p1[s], p0[s])) for s in p0 if p1[s] != p0[s]}
    return dq, dp, (Y, ij), h
st = Counter(); n = 0
for key, b0 in census():
    for b, DV in states(b0, cap=25):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind in ('BC', 'BB'):
                    if kind == 'BB' and x > y: continue
                    b1 = copy.deepcopy(b); X1 = disp(b1, x, DV[kind[0]]); Y1 = disp(b1, y, DV[kind[1]])
                    b2 = copy.deepcopy(b); Y2 = disp(b2, y, DV[kind[1]]); X2 = disp(b2, x, DV[kind[0]])
                    comm = tiles_of(b1) == tiles_of(b2)
                    inter = [s[0] for s in X1[3]] != [s[0] for s in X2[3]] or [s[1] for s in Y1[3]] != [s[1] for s in Y2[3]]
                    k = (kind, 'interacting' if inter else 'quiet')
                    st[k+('n',)] += 1; st[k+('commute',)] += comm
                    st[k+('x moved-set same',)] += (set(X1[0]), set(X1[1])) == (set(X2[0]), set(X2[1]))
                    st[k+('y moved-set same',)] += (set(Y1[0]), set(Y1[1])) == (set(Y2[0]), set(Y2[1]))
                    st[k+('x disp same',)] += X1[:2] == X2[:2]
                    st[k+('y disp same',)] += Y1[:2] == Y2[:2]
                    # union of moved sets
                    st[k+('union moved sets same',)] += (set(X1[0])|set(Y1[0]), set(X1[1])|set(Y1[1])) == (set(X2[0])|set(Y2[0]), set(X2[1])|set(Y2[1]))
                    st[k+('x,y moved 0-sets disjoint (order1)',)] += not (set(X1[0]) & set(Y1[0]))
    n += 1
    if n >= 60: break
for k in sorted(st): print(k, st[k])
