"""Step 0. The brief's control "reverse the (ii) B-list and apply iota" IS Order Reversal (c_k = iota(b_{m+1-k})), so it agrees
with (iii) exactly when G1 holds.  To decide 'bookkeeping or not', measure whether (iii) is a time-reversal of (ii) at the level
of journeys: compare the hexagon step-strings of (ii)'s k-th journey with (iii)'s (m+1-k)-th, and the A-cell each one removes."""
from common import *
from collections import Counter
st = Counter()
def journeys(b, sgn, dirn):
    u, w = b.nest[0]['u'], b.nest[0]['w']; E_ = tuple(sgn*a for a in u); N_ = tuple(sgn*a for a in w)
    fl = unique_flock(b, 0, E_, N_)
    dvec = {'E': E_, 'N': N_, 'W': tuple(-a for a in E_), 'S': tuple(-a for a in N_)}[dirn]
    key = {t: (v, box_of(b, 0, b.tiles[t][1], E_, N_)[0]) for t, v in fl.items()}
    order = sorted(fl, key=lambda t: key[t], reverse=(dirn in 'NE'))
    inv = {t: ij for ij, t in b.nest[0]['pos'].items()}; out = []
    for t in order:
        cell = inv[t]; del b.nest[0]['pos'][cell]; R, X, ij, s = b.journey(t, 0, unit_dir(dvec))
        out.append((cell, (s['a'], s['b'])))
    return out
for key, b0 in census():
    J2 = journeys(copy.deepcopy(b0), 1, 'N'); J3 = journeys(copy.deepcopy(b0), -1, 'W'); m = len(J2)
    st['cells'] += 1; st['journeys'] += m
    st['step-strings match under k<->m+1-k'] += sum(J2[k][1] == J3[m-1-k][1] for k in range(m))
    st['A-cell removed matches under k<->m+1-k'] += sum(J2[k][0] == J3[m-1-k][0] for k in range(m))
    st['A-cell removed matches under k<->k'] += sum(J2[k][0] == J3[k][0] for k in range(m))
    st['cells with all step-strings matched'] += all(J2[k][1] == J3[m-1-k][1] for k in range(m))
    st['total hexagon moves (ii)'] += sum(sum(j[1]) for j in J2); st['total hexagon moves (iii)'] += sum(sum(j[1]) for j in J3)
print(dict(st))
