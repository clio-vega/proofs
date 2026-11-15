import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585')
from rec import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
col={'tri':'#c5c5f0','sq':'#f08a8a','rh':'#48c048'}
def panel(ax, b, title, hl=()):
    for t,(k,vs) in b.tiles.items():
        z=[G.phi(P) for P in vs]; ax.fill([w.real for w in z],[w.imag for w in z],color=('#006000' if t in hl else col[k]),ec='k',lw=.6)
    for X,info in enumerate(b.nest):
        z=G.phi(info['X']); ax.text(z.real,z.imag,'ABC'[X],color='b',fontsize=14)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title)
want = ((2, 1), (2, 1), (3, 2, 1, 0), ((1, 2), (3,)))
for key, b0 in census():
    if key == want: break
x, y = (1,0), (0,1); A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
fig, axs = plt.subplots(2,2, figsize=(12,12))
panel(axs[0][0], b0, 'M', hl=(A['pos'][x], A['pos'][y]))
b1 = copy.deepcopy(b0); mig1(b1,0,x,DV['B']); panel(axs[0][1], b1, 'x->B')
b2 = copy.deepcopy(b0); mig1(b2,0,y,DV['C']); panel(axs[1][0], b2, 'y->C')
b3 = copy.deepcopy(b1); mig1(b3,0,y,DV['C']); panel(axs[1][1], b3, 'x->B ; y->C')
fig.savefig('/tmp/diamond0.png', dpi=60)
