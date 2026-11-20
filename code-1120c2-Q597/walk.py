"""PROVE 11-20 c2: print the projected walk of the rhombus R_k = s0 x s1 in each lattice, in (eta, a) coordinates.
Reads ../code-1120-Q597/dumpj5-*.pkl.   usage: python3 walk.py NJOURNEYS"""
import sys, pickle, glob
from collections import Counter
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
def rh(S): return next(t for t in S if t[0] == 'rh')
def proj(t, i): return sorted({P[2*i:2*i+2] for P in t[1]})
def show(TT, J):
    eta = ETA[TT]
    out = []
    for ty, g, n in J:
        R = rh(g); s0 = proj(R, 0); s1 = proj(R, 1)
        out.append('%s s0=%s s1=%s' % (ty, [(eta(p), p[0]) for p in s0], [(-eta(p), p[0]) for p in s1]))
    R = rh(J[-1][2]); out.append('end s0=%s s1=%s' % ([(eta(p), p[0]) for p in proj(R,0)], [(-eta(p), p[0]) for p in proj(R,1)]))
    return out
N = int(sys.argv[1]); c = 0
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for TT, M0, J in D['BC']:
        if len(J) < 6: continue
        print('T=', TT, ''.join(s[0] for s in J)); print('\n'.join(show(TT, J))); c += 1
        if c >= N: sys.exit()
