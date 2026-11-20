"""PROVE 11-20 c2: Lemma D candidates. At every step k of every dumped journey, for every removed non-rh tile t and every
unit vector w of Lambda_0 or Lambda_1 (12), classify w by dpsi = psi(t+w)-psi(t) and record whether t+w is a board tile
AFTER step k, and whether w is +-an edge vector of R_k.  D1: w edge of R_k & descending => not a tile.  D2: any descending w.
D3: split by allowed class (eta_T(w) != 0).  usage: python3 lemmaD3.py"""
import pickle, glob
from collections import Counter
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
W = [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
def rh(S): return next(t for t in S if t[0] == 'rh')
def psi(t, TT):
    eta = ETA[TT]
    return (min(-eta(P[2:]) for P in t[1]), min(eta(P[:2]) for P in t[1]))
def tr(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
def edgevecs(R):
    V = list(R[1]); return {tuple(b-a for a, b in zip(P, Q)) for P in V for Q in V}
st = Counter()
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            Bd = set(M0)
            for k, (ty, g, n) in enumerate(J):
                Bd = (Bd - set(g)) | set(n)
                EV = edgevecs(rh(g))
                for t in g:
                    if t[0] == 'rh': continue
                    pt = psi(t, TT)
                    for w in W:
                        c = 1 if any(w[:2]) else 0
                        d = psi(tr(t, w), TT); dp = (d[0]-pt[0], d[1]-pt[1])
                        e = [0, 0]; e[c] = 1
                        lab = 'asc' if dp == tuple(e) else ('desc' if dp == tuple(-x for x in e) else 'flat')
                        wl = w[:2] if any(w[:2]) else w[2:]
                        allowed = ETA[TT](wl) != 0
                        st[(fam, lab, 'w edge of R_k', w in EV, 'allowed class', allowed, 't+w tile after k', tr(t, w) in Bd)] += 1
for k in sorted(st, key=str): print(st[k], k)
