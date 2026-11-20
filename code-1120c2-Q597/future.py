"""PROVE 11-20 c2: universal test of Theorem G1c-future on every dumped journey (B/C and B/B-second).
For every step j, every removed non-rh tile u of H_j, every edge vector v of R_j (both signs, both lattices):
  'ascending'  : psi(u+v) = psi(u) + e_i, i = the lattice of v  -> claim: u+v is removed by NO step k > j.
  'descending' : psi(u+v) = psi(u) - e_i  (lesion; no claim).
Also Lemma W: if a removed tile has psi_i < mu_i then psi_i = mu_i - 1 and mu_{j+1,i} = mu_{j,i} - 1.
usage: python3 future.py"""
import pickle, glob
from collections import Counter
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
def rh(S): return next(t for t in S if t[0] == 'rh')
def psi(t, TT):
    eta = ETA[TT]
    return (min(-eta(P[2:]) for P in t[1]), min(eta(P[:2]) for P in t[1]))
def tr(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
def edgevecs(R):
    V = list(R[1]); out = set()
    for P in V:
        for Q in V:
            d = tuple(b-a for a, b in zip(P, Q))
            nz0, nz1 = any(d[:2]), any(d[2:])
            if (nz0 != nz1) and max(abs(x) for x in d) <= 1: out.add(d)   # unit vector in one lattice
    return out
st = Counter()
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            mu = [psi(rh(g), TT) for _, g, _ in J] + [psi(rh(J[-1][2]), TT)]
            remat = {}
            for k, (_, g, _) in enumerate(J):
                for t in g:
                    if t[0] != 'rh': remat.setdefault(t, []).append(k)
            for j, (ty, g, n) in enumerate(J):
                R = rh(g); EV = edgevecs(R)
                st[(fam, 'edge vectors of R per step', len(EV))] += 1
                for u in g:
                    if u[0] == 'rh': continue
                    pu = psi(u, TT)
                    for i in (0, 1):
                        if pu[i] < mu[j][i]:
                            st[('LemmaW', fam, 'psi_i = mu_i - 1', pu[i] == mu[j][i]-1, 'mu drops by 1 at j', mu[j+1][i] == mu[j][i]-1)] += 1
                    for v in EV:
                        i = 0 if any(v[:2]) else 1     # v in Lambda_0 changes psi_1; v in Lambda_1 changes psi_0
                        comp = 1 if i == 0 else 0
                        w = tr(u, v); pw = psi(w, TT)
                        dpsi = (pw[0]-pu[0], pw[1]-pu[1])
                        e = [0, 0]; e[comp] = 1
                        lab = 'ascending' if dpsi == tuple(e) else ('descending' if dpsi == tuple(-x for x in e) else 'other %s' % (dpsi,))
                        later = [k for k in remat.get(w, []) if k > j]
                        st[(fam, TT, lab, 'u+v removed by a later step', bool(later))] += 1
                        if lab == 'ascending':
                            st[(fam, 'ascending: psi_comp(u+v) - mu_{j+1}', pw[comp] - mu[j+1][comp])] += 1
for k in sorted(st, key=str): print(st[k], k)
