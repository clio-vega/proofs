"""Draw M_P and Mx_P (and successive boards along both middles) near the middle, for one record per catalogue config."""
import sys, pickle, cmath, math
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
from mosaic import Geo
from replay import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
G = Geo(-1)
COL = {'rh': '#d62728', 'sq': '#9ecae1', 'tri': '#fdd0a2'}
def order(vs):
    pts = [G.phi(P) for P in vs]; c = sum(pts)/len(pts)
    return sorted(pts, key=lambda z: cmath.phase(z-c))
def draw(ax, Bd, hl, title, reg):
    for k, vs in Bd:
        pts = order(vs)
        if not any(abs(z-reg) < 3.2 for z in pts): continue
        t = (k, vs)
        ax.add_patch(Polygon([(z.real, z.imag) for z in pts], closed=True, fc=COL[k] if k != 'tri' else ('#fdd0a2' if len({v[2:] for v in vs}) == 1 else '#c7e9c0'),
                             ec='k', lw=2.5 if t in hl else 0.5, alpha=0.95))
    ax.set_xlim(reg.real-3, reg.real+3); ax.set_ylim(reg.imag-3, reg.imag+3); ax.set_aspect('equal'); ax.set_title(title, fontsize=7); ax.axis('off')
want = {(0, 1): 2, (1, 1): 4}; got = {}
for r in R:
    if r['kind'] != 'BC' or not r['Us'] or len(r['Us']) != 1: continue
    side, A, B, Jo = clean(r); nU = r['Us'][0][1]
    if (side, nU) not in want: continue
    mA, mB, p, q = middles(A, B)
    o = tuple(-v for v in cen(rh(mA[0][1])))
    K = (tuple(tr(s, o) for s in mA), tuple(tr(s, o) for s in mB))
    got.setdefault((side, nU), {}).setdefault(K, r)
for (side, nU), d in got.items():
    for i, (K, r) in enumerate(d.items()):
        _, A, B, Jo = clean(r); mA, mB, p, q = middles(A, B)
        M0 = board0(r); BA = apply(M0, A[:p]); BB = apply(apply(M0, Jo), B[:p])
        reg = G.phi(cen(rh(mA[0][1])))
        n = len(mA); fig, axs = plt.subplots(2, n+1, figsize=(2.4*(n+1), 5))
        for j in range(n+1):
            D = BA ^ BB
            draw(axs[0, j], BA, D, f'M step{j} ' + (mA[j][0] if j < n else 'end'), reg)
            draw(axs[1, j], BB, D, f'Mx step{j} ' + (mB[j][0] if j < n else 'end'), reg)
            if j < n: BA = apply(BA, [mA[j]]); BB = apply(BB, [mB[j]])
        fig.suptitle(f'side {side} |U|={nU} config {i} key {r["key"]}  (thick = D)', fontsize=8)
        fig.savefig(f'viz-{side}-{nU}-{i}.png', dpi=110); plt.close(fig)
        print(side, nU, i, [s[0] for s in mA], [s[0] for s in mB])
