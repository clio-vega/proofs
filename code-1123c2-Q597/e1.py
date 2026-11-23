"""Retreat (case b): (E1) J_o's step e+1 is not (1,a_T').  Mechanism test: does the (1,a_T') T'-configuration at J_o's image
R^o_{e+1} of the split contain a non-rh tile created by J_o at e (then no-revisit forbids it)?  Same for the other types
(control: the types that DO occur must not meet Cre_J(e)).  Census: J_o's step-e+1 types in case b.   usage: python3 e1.py"""
import sys, os, pickle
from collections import Counter
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
nr = lambda S: {u for u in S if u[0] != 'rh'}
for T in 'BC':
    Tp = OTHER[T]; Rs = CF[T][(k0, 'b')][0]
    HJ, RJe1, CJ = step(Tp, Rs, 'b')
    for t1 in ('a0', 'a1', 'b'):
        H1, _, _ = step(Tp, RJe1, t1)
        print('T', T, 'J_o (toward', Tp + ') step e+1 type (1,%s)' % t1, 'meets Cre_J(e) non-rh:', len(nr(H1) & nr(CJ)))
