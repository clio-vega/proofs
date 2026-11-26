"""Provenance of the tiles removed by the front-window steps E_1 (on M0), E^L_0, E^L_1 (on M_L), L_1 (on M0).
Labels: 'L<j>' created by J_L's step p_L+j (relative); 'EL<j>' created by E^L's step p_E+j; 'E<j>' created by J_E's step;
'M0:inH(E0)|inH(L0)|other' an M0 tile, tagged by whether it lies in a contact hexagon."""
import pickle, glob
from collections import Counter
st = Counter()
def trace(runs):
    """runs: list of (label_prefix, J, p). returns dict tile -> label of the step creating it (latest)."""
    made = {}
    for lab, J, p in runs:
        for i, s in enumerate(J):
            for t in s[1]: made.pop(t, None)
            for t in s[2]: made[t] = '%s%+d' % (lab, i - p)
    return made
for fn in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(fn, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        HE0, HL0 = set(AE[pE][1]), set(AL[pL][1])
        def lab(t, made):
            if t in made: return made[t]
            return 'M0:' + ('H(E0)&H(L0)' if t in HE0 and t in HL0 else 'H(E0)' if t in HE0 else 'H(L0)' if t in HL0 else 'other')
        # E^L steps: board = M0 . J_L . E^L[:j]
        for j in (0, 1):
            made = trace([('L', AL, pL), ('EL', ArE[:pE+j], pE)])
            st[(case, 'E^L_%d' % j, tuple(sorted(lab(t, made) for t in ArE[pE+j][1])))] += 1
        made = trace([('E', AE[:pE+1], pE)])
        st[(case, 'E_1', tuple(sorted(lab(t, made) for t in AE[pE+1][1])))] += 1
        made = trace([('L', AL[:pL+1], pL)])
        st[(case, 'L_1', tuple(sorted(lab(t, made) for t in AL[pL+1][1])))] += 1
for k in sorted(st, key=str): print(st[k], k)
