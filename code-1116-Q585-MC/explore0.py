"""Exploration: decompose tiles by Z^4 parts. 0-triangles (all 1-parts equal), 1-triangles, parallelograms.
Q: are 0-parts of 0-triangles a fixed set (independent of state)?  Exact integer arithmetic only."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1115c2-Q585-G1')
from rec import *
from collections import Counter
def decomp(b):
    T0, T1, P = Counter(), Counter(), Counter()
    for kind, vs in b.tiles.values():
        p0 = {v[:2] for v in vs}; p1 = {v[2:] for v in vs}
        if kind == 'tri':
            if len(p1) == 1: T0[frozenset(p0)] += 1
            elif len(p0) == 1: T1[frozenset(p1)] += 1
            else: raise ValueError('mixed triangle')
        else:
            assert len(p0) == 2 and len(p1) == 2; P[(frozenset(p0), frozenset(p1))] += 1
    return T0, T1, P
n = 0; st = Counter()
for key, b0 in census():
    T0a, T1a, Pa = decomp(b0)
    st['multiplicity>1 in T0'] += any(v > 1 for v in T0a.values())
    for b, DV in states(b0, cap=15):
        T0, T1, Pz = decomp(b)
        st['states'] += 1; st['T0 set same as start'] += T0 == T0a; st['T1 same'] += T1 == T1a; st['P same'] += Pz == Pa
    n += 1
    if n >= 30: break
print(st)
