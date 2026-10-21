import sys; sys.path.insert(0,'/home/clio/projects/scratch/20261021c2')
from g3 import *

print('=== MEASUREMENT A: ord_S from the SYMMETRISATION definition (independent of the')
print('    root product entirely).  Generic point: q and w_{j!=m} symbolic.')
for N in (1,2,3):
    w=wvars(N)
    for X in ('B','C'):
        for eps in (1,-1):
            agree=tot=0; bad=[]
            for I in subsets(N):
                A=Atilde(N,w,I,I,X=X)
                for m in range(1,N+1):
                    o=ord_S(A,w,m,eps); pred=1 if m in I else 0
                    tot+=1
                    if o==pred: agree+=1
                    else: bad.append((I,m,o,pred))
            print(f'  N={N} X={X} eps={eps:+d}: {agree:3d}/{tot:3d} = [m in I]', flush=True)
            if bad:
                from collections import Counter
                print('      disagreements (I,m,measured,predicted):', bad[:6],
                      ' pattern:', dict(Counter((o,p) for _,_,o,p in bad)), flush=True)
