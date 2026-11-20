"""PROVE 11-20 c2: per-step changes of (min eta, min a) of the rhombus projections s0, s1, on every dumped journey.
usage: python3 latstat.py"""
import pickle, glob
from collections import Counter
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
def rh(S): return next(t for t in S if t[0] == 'rh')
def proj(t, i): return {P[2*i:2*i+2] for P in t[1]}
st = Counter()
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            eta = ETA[TT]; mono = [True]*4
            for ty, g, n in J:
                R, Rn = rh(g), rh(n)
                for i in (0, 1):
                    s, sn = proj(R, i), proj(Rn, i)
                    de = min(eta(p) for p in sn) - min(eta(p) for p in s)
                    da = min(p[0] for p in sn) - min(p[0] for p in s)
                    st[(fam, TT, 'lattice', i, 'type', ty, 'd min eta', de, 'd min a', da)] += 1
                    if da < 0: mono[2*i] = False
                    if da > 0: mono[2*i+1] = False
            st[(fam, TT, 'journey: min a of s0 nondecr, nonincr, of s1 nondecr, nonincr', tuple(mono))] += 1
for k in sorted(st, key=str): print(st[k], k)
