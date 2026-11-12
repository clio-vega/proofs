from flock import *
import time
lam=(2,1); k0,k1=3,3
for nu in [(1,1),(2,)]:
    mu=(1,)
    boxes={(i,j) for j,l in enumerate(lam) for i in range(l)}-{(0,0)}
    for f in lr_tableaux(boxes,nu):
        t=time.time()
        g,info,b=cor34(f,lam,k0,k1,verbose=True)
        print('f',f,'-> G\'',g,'shape',skew_of(set(g)),'content',content(g),'F\'',info[2],skew_of(set(info[2])),'%.1fs'%(time.time()-t))
