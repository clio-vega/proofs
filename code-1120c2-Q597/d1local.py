"""PROVE 11-20 c2: is Lemma D1 local?  For every step k, removed non-rh t, descending edge vector w of R_k:
(i) centroid of phi(t+w) strictly inside the hexagon H_k;  (ii) t+w is not a created tile of step k.
If (i) and (ii) hold, t+w overlaps the new tiling of H_k in its interior and so is not a board tile after step k.
Also evaluated on the abstract good configurations of localenum (board-free).   usage: python3 d1local.py"""
import sys, pickle, glob, cmath, math
from collections import Counter
sys.path.insert(0, '../code-1120-Q597')
from psilib import psi, rh, tr, edgevecs, ETA
ZETA = cmath.exp(1j*math.pi/3); z0 = cmath.exp(1j*math.pi/12); z1 = cmath.exp(-1j*math.pi/12)
def phi(P): return z0*(float(P[0]) + float(P[1])*ZETA) + z1*(float(P[2]) + float(P[3])*ZETA)
def hull(pts):
    pts = sorted(set((round(p.real, 9), round(p.imag, 9)) for p in pts))
    def cr(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 1e-9: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 1e-9: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]
def strictly_inside(z, H):
    n = len(H)
    for i in range(n):
        a, b = H[i], H[(i+1) % n]
        if (b[0]-a[0])*(z.imag-a[1]) - (b[1]-a[1])*(z.real-a[0]) <= 1e-9: return False
    return True
def test(TT, g, n, st, fam):
    R = rh(g); H = hull([phi(P) for t in g for P in t[1]])
    for t in g:
        if t[0] == 'rh': continue
        pt = psi(t, TT)
        for w in edgevecs(R):
            if not any(w): continue
            c = 1 if any(w[:2]) else 0
            s = tr(t, w); ps = psi(s, TT)
            if ps[c] - pt[c] != -1: continue
            cz = sum(phi(P) for P in s[1]) / len(s[1])
            st[(fam, 'desc R-edge translate: centroid strictly inside H', strictly_inside(cz, H), 'is a created tile', s in n)] += 1
if __name__ == '__main__':
    st = Counter()
    for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
        D = pickle.load(open(fn, 'rb'))
        for fam in ('BC', 'BB2'):
            for TT, M0, J in D[fam]:
                for ty, g, n in J: test(TT, g, n, st, fam)
    for k in sorted(st, key=str): print(st[k], k)
