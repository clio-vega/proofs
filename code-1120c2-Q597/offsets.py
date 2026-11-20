"""PROVE 11-20 c2: the 4-offset of every removed non-rh tile relative to the rhombus, on every dumped journey.
Om(t, R) = (psi_0(t)-mu_0, psi_1(t)-mu_1, A_0(t)-lam_0, A_1(t)-lam_1), A_i = min first lattice coordinate over pi_i,
lam_i = A_i(R).  Prints the set Q of offsets of removed tiles, per (T, tile kind).   usage: python3 offsets.py"""
import pickle, glob
from collections import Counter
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
def rh(S): return next(t for t in S if t[0] == 'rh')
def kind(t):
    if t[0] != 'tri': return t[0]
    return 'tri0' if len({P[2:] for P in t[1]}) == 1 else 'tri1'
def vec4(t, TT):
    eta = ETA[TT]
    return (min(-eta(P[2:]) for P in t[1]), min(eta(P[:2]) for P in t[1]), min(P[0] for P in t[1]), min(P[2] for P in t[1]))
def om(t, R, TT): a, b = vec4(t, TT), vec4(R, TT); return tuple(int(x - y) for x, y in zip(a, b))
if __name__ == '__main__':
    st = Counter()
    for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
        D = pickle.load(open(fn, 'rb'))
        for fam in ('BC', 'BB2'):
            for TT, M0, J in D[fam]:
                for ty, g, n in J:
                    R = rh(g)
                    for t in g:
                        if t[0] == 'rh': continue
                        st[(TT, kind(t), om(t, R, TT))] += 1
    for k in sorted(st, key=str): print(st[k], k)
