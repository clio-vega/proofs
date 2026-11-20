"""PROVE 11-20 c2: live control for Theorem F.  Same as future2.py, but w ranges over ALL 12 unit vectors.
Theorem F claims: w an edge vector of R_j and ascending => t+w never removed after j.  The control asks whether the
claim FAILS when w is ascending but NOT an edge vector of R_j (it should, else the hypothesis is idle).
Labels asc/flat/desc by psi_c(u+w)-psi_c(u); flat and desc are the controls.  usage: python3 fcontrol.py"""
import pickle, glob
from collections import Counter
from psilib import psi, rh, tr, edgevecs, ETA
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
W = [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
st = Counter()
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            remat = {}
            for k, (_, g, _) in enumerate(J):
                for t in g:
                    if t[0] != 'rh': remat.setdefault(t, []).append(k)
            for j, (ty, g, n) in enumerate(J):
                EV = edgevecs(rh(g))
                for u in g:
                    if u[0] == 'rh': continue
                    pu = psi(u, TT)
                    for w in W:
                        c = 1 if any(w[:2]) else 0
                        t = tr(u, w); pt = psi(t, TT)
                        lab = {1: "asc", 0: "flat", -1: "desc"}.get(pt[c] - pu[c], "other")
                        later = any(k > j for k in remat.get(t, []))
                        st[(fam, lab, "edge of R_j", w in EV, "u+w removed later", later)] += 1
for k in sorted(st, key=str): print(st[k], k)
