"""Gashed-puzzle frames (psi-images) of the four journeys of the N-th (SW) instance.
Each step frame: triangles/squares mapped by psi=(P0+P2,P1+P3) in the triangular lattice; 0-edges blue, 1-edges red,
square diagonal (10) dotted; moving rhombus = thick black segment (the gash).  usage: python3 gashpic.py N"""
import sys
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = complex(0.5, 3**.5/2)
def z(P): return (P[0]+P[2]) + (P[1]+P[3])*W
FR = []
_s2 = migrate.Board.step
def snap(self, R, d):
    R2 = _s2(self, R, d)
    if R2 is not None and REC[0]: FR.append((copy.deepcopy(self.tiles), R2, REC[0]))
    return R2
migrate.Board.step = snap; REC = ['']
def draw(ax, tiles, R, title):
    for t, (k, vs) in tiles.items():
        n = len(vs)
        for i in range(n):
            P, Q = vs[i], vs[(i+1) % n]; a, b = z(P), z(Q)
            d = [q-p for p, q in zip(P, Q)]
            if k == 'rh':
                if t == R: ax.plot([a.real, b.real], [a.imag, b.imag], 'k-', lw=3.5)
                continue
            c = 'b' if not (d[2] or d[3]) else ('r' if not (d[0] or d[1]) else 'g')
            ax.plot([a.real, b.real], [a.imag, b.imag], c+'-', lw=1)
        if k == 'sq':
            # short diagonal: the two vertices whose difference is mixed
            pts = list(vs)
            for i in range(4):
                for j in range(i+1, 4):
                    d = [q-p for p, q in zip(pts[i], pts[j])]
                    if (d[0] or d[1]) and (d[2] or d[3]) and abs(z(pts[i])-z(pts[j])) < 1.01:
                        ax.plot([z(pts[i]).real, z(pts[j]).real], [z(pts[i]).imag, z(pts[j]).imag], 'g:', lw=1)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title, fontsize=7)
N = int(sys.argv[1]); cnt = 0
for key, b0 in census():
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            if cnt < N: cnt += 1; continue
            FR.clear()
            init = copy.deepcopy(b.tiles)
            REC[0] = 'K x->B'; b1 = copy.deepcopy(b); mv(b1, x, dirs['B']); REC[0] = 'G y->C'; mv(b1, y, dirs['C'])
            REC[0] = 'H x->C'; b2 = copy.deepcopy(b); mv(b2, x, dirs['C']); REC[0] = 'Gp y->B'; mv(b2, y, dirs['B'])
            REC[0] = ''
            print(key, x, y, len(FR), 'same', tiles_of(b1) == tiles_of(b2))
            byj = {}
            for f in FR: byj.setdefault(f[2], []).append(f)
            ncol = max(len(v) for v in byj.values()) + 1
            fig, axs = plt.subplots(4, ncol, figsize=(2.6*ncol, 10.4))
            for r, nm in enumerate(['K x->B', 'G y->C', 'H x->C', 'Gp y->B']):
                for c in range(ncol): axs[r][c].axis('off')
                if r in (0, 2): draw(axs[r][0], init, None, 'tau')
                for c, (tl, R, _) in enumerate(byj[nm]): draw(axs[r][c+1], tl, R, f'{nm} step {c+1}')
            fig.savefig(f'/home/clio/projects/proofs/code-1128c2-Q585-SW/gash-{N}.png', dpi=55, bbox_inches='tight'); sys.exit()
