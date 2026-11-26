"""(Merge) J_E, J_L share identical steps on M0; (Merge') E^L, L^E share identical steps (on different boards M_L, M_E).
Relative indices to contacts; corridor length k; and is the re-run corridor = journey corridor translated by one vector?"""
import pickle, glob, sys
from collections import Counter
from match import canon, vsub
st = Counter()
def K(s): return (s[0], s[1], s[2])
def shared(A, B, pa, pb):
    idx = {K(s): i for i, s in enumerate(B)}
    return [(i - pa, idx[K(s)] - pb) for i, s in enumerate(A) if K(s) in idx]
for f in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(f, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        s0 = shared(AE, AL, pE, pL); s1 = shared(ArE, ArL, pE, pL)
        k = laL - 3 if case == 'b' else laL - 2
        st[(case, 'k', k, 'M0 shared', tuple(s0), 'rerun shared', tuple(s1))] += 1
        # translation of corridor between journey and re-run (E's corridor steps vs E^L's shared steps)
        if s0 and s1:
            vs = set()
            for (i0, j0), (i1, j1) in zip(s0, s1):
                a, ma = canon(AE[pE+i0]); b, mb = canon(ArE[pE+i1])
                vs.add(vsub(mb, ma) if a == b else None)
            st[(case, 'T_L', TL, 'corridor re-run = corridor + one vector', tuple(vs))] += 1
for k in sorted(st, key=str): print(st[k], k)
