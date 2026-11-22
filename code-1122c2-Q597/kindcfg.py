"""Local trichotomy on the configuration table: J_o toward T' (other target), alpha = delta_T (T = clean target).
For every good T'-configuration (class, type) and every removed non-rh tile u: is u - alpha CREATED by the same step
(kind B, local)?  is -alpha an edge vector of u (A-shaped)?  kind of tile.  usage: python3 kindcfg.py"""
import sys, os
from fractions import Fraction as Fr
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
exec(open('/home/clio/projects/proofs/code-1122c2-Q597/sidecfg.py').read().split('st = Counter()')[0])
A_CLASS = frozenset({(0, 1), (0, -1)})
DELTA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
def edge_of(t, v):
    V = list(t[1]); return any(tuple(Fr(b)-Fr(a) for a, b in zip(P, Q)) == tuple(Fr(x) for x in v) for P in V for Q in V)
for T, Tp in (('B', 'C'), ('C', 'B')):
    ma = tuple(-x for x in DELTA[T])
    C = Wd.configs(Tp)
    for (ck, typ), (R, rem, R2, cre) in sorted(C.items(), key=str):
        cl = 0 if ck == A_CLASS else 1
        for u in sorted(rem, key=str):
            if u[0] == 'rh': continue
            um = shift(u, ma)
            print('clean', T, 'Jo', Tp, 'class', cl, typ, 'tile', u[0], 'u-alpha created by step (B)', um in cre,
                  'A-shaped', edge_of(u, ma), 'u-alpha in H (removed)', um in rem)
