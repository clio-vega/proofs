import sys; sys.path.insert(0,'.')
from tm import *
from fibre import gvars
import sympy as sp
def sig_tm(h,w,nmax):
    g=gvars(h,w); S,A=transfer(h,w,g); idx={m:a for a,m in enumerate(S)}
    u=sp.zeros(len(S),1); u[idx[tuple([0]*(h-1)+[w])],0]=1
    e=sp.Matrix([[eps(m) for m in S]])
    out=[]; v=u
    for n in range(1,nmax+1):
        v=sp.expand(A*v); out.append(sp.expand((e*v)[0,0]))
    return out
h,w=4,2
S=sig_tm(h,w,14)
b=[S[2*r-1] for r in range(1,8)]   # Sigma_2,4,...,14
for sz in [2,3,4]:
    H=sp.Matrix(sz,sz, lambda i,j: b[i+j])
    d=sp.expand(H.det())
    print(f'h={h} w={w}: Hankel {sz}x{sz} on (Sigma_2..Sigma_{2*(2*sz-1)}): det==0 ? {d==0}', flush=True)
