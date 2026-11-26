"""Identical steps (same tile sets) shared by J_E and J_L on M0: the corridor.  Also E's re-run vs L's re-run."""
import pickle, glob, sys
from collections import Counter
sys.argv=['x','5','1','24']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
st = Counter(); shown = 0
def key(s): return (s[0], s[1], s[2])
for f in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(f, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        idxL = {key(s): i for i, s in enumerate(AL)}
        shared = [(i, idxL[key(s)]) for i, s in enumerate(AE) if key(s) in idxL]
        offs = sorted(set(j - i for i, j in shared))
        rel = [(i - pE, j - pL) for i, j in shared]
        st[(case, 'shared steps rel to contacts (E,L)', tuple(rel))] += 1
        if shown < 6:
            shown += 1
            print(case, 'E->', TE, 'pE', pE, ws(typed(AE)), '|| L->', TL, 'pL', pL, ws(typed(AL)), 'shared', shared)
for k in sorted(st, key=lambda k: -st[k]): print(st[k], k)
