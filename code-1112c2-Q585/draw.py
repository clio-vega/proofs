import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from flock import *
def draw(b,fn,rot=1,labels=None):
    fig,ax=plt.subplots(figsize=(7,7)); col={'tri':'#c5c5f0','sq':'#f08a8a','rh':'#48c048'}
    for t,(k,vs) in b.tiles.items():
        z=[G.phi(P)*rot for P in vs]; ax.fill([w.real for w in z],[w.imag for w in z],color=col[k],ec='k',lw=.6)
        if labels and t in labels:
            c=sum(z)/len(z); ax.text(c.real,c.imag,str(labels[t]),ha='center',va='center')
    for X,info in enumerate(b.nest):
        z=G.phi(info['X'])*rot; ax.text(z.real,z.imag,'ABC'[X],color='b',fontsize=14)
    ax.set_aspect('equal'); ax.axis('off'); fig.savefig(fn,dpi=70)
if __name__=='__main__':
    # Fig 7 analogue: n=8,d=5, beta=(2,2,1)
    b=canonical((2,2,1),3,5); draw(b,'/tmp/fig7a.png')
    b=canonical((2,2,1),5,3); draw(b,'/tmp/fig7b.png')
