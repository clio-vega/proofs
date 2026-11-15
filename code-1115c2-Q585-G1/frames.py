"""Frame-by-frame trace of y's C-journey (solo, and after x->B) in the smallest interacting diamond.
Highlight: moving rhombus (dark green), the hexagon being rotated (thick outline), tiles of x's B-wake (yellow)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585')
from rec import *
import migrate
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
col={'tri':'#c5c5f0','sq':'#f08a8a','rh':'#48c048'}
want = ((2, 1), (2, 1), (3, 2, 1, 0), ((1, 2), (3,)))
for key, b0 in census():
    if key == want: break
A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}; x, y = (1,0), (0,1)
def frames(b, ij, dvec, wake, row, axs, label):
    t = b.nest[0]['pos'].pop(ij); d = unit_dir(dvec); R = t; k = 0
    while True:
        vs = b.tiles[R][1]; done = False
        for X, info in enumerate(b.nest):
            if X == 0: continue
            p = b.packed_at(info, vs)
            if p is not None and p not in info['pos'] and b.is_ideal(set(info['pos']) | {p}): done = True
        ax = axs[row][k]
        before = dict(b.tiles)
        R2 = None if done else b.step(R, d)
        gone = set(before) - set(b.tiles) if R2 is not None else set()
        for tt, (kk, vv) in before.items():
            z = [G.phi(P) for P in vv]
            c = '#006000' if tt == R else ('#f0e060' if frozenset(vv) in wake else col[kk])
            ax.fill([w.real for w in z], [w.imag for w in z], color=c, ec='k', lw=(2.2 if tt in gone else .4))
        ax.set_aspect('equal'); ax.axis('off'); ax.set_title(f'{label} {k}' + (f' ({b.last})' if R2 is not None else ' landed'), fontsize=8)
        pass
        if done: return k
        R = R2; k += 1
b1 = copy.deepcopy(b0); pre = {frozenset(v) for _, v in b1.tiles.values()}
mig1(b1, 0, x, DV['B']); post = {frozenset(v) for _, v in b1.tiles.values()}
wake = post - pre
fig, axs = plt.subplots(2, 7, figsize=(28, 9))
for a in axs.flat: a.axis('off')
frames(copy.deepcopy(b0), y, DV['C'], set(), 0, axs, 'solo')
frames(copy.deepcopy(b1), y, DV['C'], wake, 1, axs, 'after x->B')
fig.tight_layout(); fig.savefig('/tmp/frames.png', dpi=45)
print(len(wake), 'wake tiles')
