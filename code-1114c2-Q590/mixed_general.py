"""General mixed commutation: from EVERY state reachable by mixed single migrations (any corner, either direction),
for every ordered pair of DISTINCT current corners (x to B, y to C): x->B then y->C  vs  y->C then x->B.
Non-triviality statistic: diamonds where some journey's hexagon step-string differs between the two paths (the moves interact).
Control (same instrument): same-direction pairs (x->B, y->B)."""
from common import *
from collections import Counter
st = Counter(); bad = []
def key_of(b): return (frozenset(b.nest[0]['pos']), tuple(tiles_of(b)))
for key, b0 in census():
    A = b0.nest[0]; uA, wA = A['u'], A['w']; DV = {'B': wA, 'C': uA}
    stack = [b0]; seen = set(); nstates = 0
    while stack and nstates < 60:
        b = stack.pop(); k = key_of(b)
        if k in seen: continue
        seen.add(k); nstates += 1; cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind in ('BC', 'BB'):
                    if kind == 'BB' and x > y: continue
                    b1 = copy.deepcopy(b); _, l1x, s1x = mig1(b1, 0, x, DV[kind[0]]); _, l1y, s1y = mig1(b1, 0, y, DV[kind[1]])
                    b2 = copy.deepcopy(b); _, l2y, s2y = mig1(b2, 0, y, DV[kind[1]]); _, l2x, s2x = mig1(b2, 0, x, DV[kind[0]])
                    st[kind+' diamonds'] += 1
                    ok = tiles_of(b1) == tiles_of(b2) and (l1x, l1y) == (l2x, l2y)
                    st[kind+' commute'] += ok
                    inter = dict(s1x) != dict(s2x) or dict(s1y) != dict(s2y)
                    st[kind+' interacting (step-strings differ)'] += inter
                    st[kind+' interacting AND commute'] += inter and ok
                    if kind == 'BC' and not ok and len(bad) < 2: bad.append((key, x, y))
    st['cells'] += 1
print(dict(st)); print(bad)
