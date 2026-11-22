"""PROVE 11-22 c2: kinds of the wake tiles read at step p+1, on the four (CC) contact configurations (placement as in
code-1121c2-Q597/hsderive.py).  For each tile of A's b-hexagon at R_{p+1} that is a tile of J's step c0 or c0+1:
B-local (u - alpha created by that J step), A-shaped (-alpha an edge vector of u), else.  usage: python3 hskind.py"""
import os, sys
from fractions import Fraction as Fr
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nfor cc in CCS:')[0]
exec(src)
DELTA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
def edge_of(t, v):
    V = list(t[1]); return any(tuple(Fr(b)-Fr(a) for a, b in zip(P, Q)) == tuple(Fr(x) for x in v) for P in V for Q in V)
for (T, t0, t1, tA, off) in CCS:
    Tp = OTHER[T]; ma = tuple(-x for x in DELTA[T])
    RJ0 = canon(Tp, t0[0]); H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1])
    RAc = canon(T, tA[0]); sh = tuple(b + o - a for a, b, o in zip(base(RAc), base(RJ0), off)); RAp = trs(RAc, sh)
    Hp, RAp1, Cp = step(T, RAp, tA[1])
    al = tuple(b - a for a, b in zip(base(RAp), base(RAp1)))
    Hb, Rb, Cb = step(T, RAp1, 'b')
    for u in sorted(Hb, key=str):
        if u[0] == 'rh': continue
        for lab, H, Cr, tj in (('c0', H0, C0, t0), ('c0+1', H1, C1, t1)):
            if u in H:
                um = trs(u, ma)
                print(T, 'contact', t0, t1, '| tile', u[0], 'from J', lab, tj, '| B-local', um in Cr, '| A-shaped', edge_of(u, ma), '| alpha check', al == tuple(Fr(x) for x in DELTA[T]))
