"""Per triangle k: which of J_L, J_E, J_L^E (L re-run), J_E^L (E re-run) move it, and the displacement each applies.
G3' per triangle:  d_L + d_{E^L} = d_E + d_{L^E}.  Pattern = (side rel. to s, tuple of displacements or '-')."""
import pickle, glob, sys
from collections import Counter
from lookstep import delta
st = Counter()
def moves(J):
    D = {}
    for s in J:
        for k, v in delta(s).items(): D[k] = v
    return D
for f in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(f, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        s = 1 if TL == 'B' else 0
        Ds = [moves(J) for J in (AL, ArE, AE, ArL)]   # L, E^L | E, L^E
        keys = set().union(*Ds)
        bad = 0
        for k in keys:
            d = [D.get(k, '-') for D in Ds]
            st[(case, 'T_L', TL, 'side', 's' if k[0] == s else 's\'', 'L,E^L | E,L^E', tuple(d))] += 1
for k in sorted(st, key=str): print(st[k], k)
