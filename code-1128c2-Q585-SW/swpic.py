"""Draw the n-th (SW) column instance on the size-5 closure: tau, tau;x->B, tau;x->C, L, R, with each journey's
created tiles shaded (wakes).  usage: python3 swpic.py N"""
import sys; sys.argv += ['0']
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
G = migrate.G
def draw(ax, b, wakes=(), title=''):
    col = {'tri': '#e8e8f8', 'sq': '#f6c8c8', 'rh': '#48c048'}
    wk = {}
    for c, J in wakes:
        for st in J:
            for vs in st: wk[frozenset(vs)] = c
    for t, (k, vs) in b.tiles.items():
        z = [G.phi(P) for P in vs]; fc = col[k]
        if k != 'rh' and frozenset(vs) in wk: fc = wk[frozenset(vs)]
        ax.fill([w.real for w in z], [w.imag for w in z], color=fc, ec='k', lw=.5)
    for X, info in enumerate(b.nest):
        z = G.phi(info['X']); ax.text(z.real, z.imag, 'ABC'[X], color='b', fontsize=12)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title, fontsize=8)
N = int(sys.argv[1]); cnt = 0
for key, b0 in census():
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            if cnt < N: cnt += 1; continue
            b1 = copy.deepcopy(b); lxB, Jx1 = mv(b1, x, dirs['B']); bxB = copy.deepcopy(b1); lyC, Jy1 = mv(b1, y, dirs['C'])
            b2 = copy.deepcopy(b); lxC, Jx2 = mv(b2, x, dirs['C']); bxC = copy.deepcopy(b2); lyB, Jy2 = mv(b2, y, dirs['B'])
            print(key, 'x', x, 'y', y, 'A-nest', sorted(S)); print(' L: x->B', lxB, len(Jx1), 'steps; y->C', lyC, len(Jy1))
            print(' R: x->C', lxC, len(Jx2), 'steps; y->B', lyB, len(Jy2)); print(' same', tiles_of(b1) == tiles_of(b2))
            print(' suffix match C-side', suffix_len(Jy1, Jx2), 'B-side', suffix_len(Jy2, Jx1))
            fig, axs = plt.subplots(1, 5, figsize=(25, 5.5))
            draw(axs[0], b, title=f'tau  x={x} y={y}')
            draw(axs[1], bxB, [('#9fd0ff', Jx1)], 'x->B (blue wake)')
            draw(axs[2], bxC, [('#ffd27f', Jx2)], 'x->C (orange wake)')
            draw(axs[3], b1, [('#9fd0ff', Jx1), ('#c9a0ff', Jy1)], 'L: x->B, y->C (purple)')
            draw(axs[4], b2, [('#ffd27f', Jx2), ('#a0e0a0', Jy2)], 'R: x->C, y->B (green)')
            fig.savefig(f'/home/clio/projects/proofs/code-1128c2-Q585-SW/swpic-{N}.png', dpi=60, bbox_inches='tight'); sys.exit()
