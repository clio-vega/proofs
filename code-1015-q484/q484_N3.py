import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp, time, q484_engine as E
from q484_struct import build, solve_c
q=E.q
N=3
t0=time.time(); w,S,A=build(N); print(f'matrix built [{time.time()-t0:.0f}s]',flush=True)
def is_laurent(f):
    f=sp.cancel(sp.together(f)); n,d=sp.fraction(f)
    return len(sp.Poly(sp.expand(d),*w).monoms())==1
nA=totA=nC=totC=nP=totP=0; vA=[]; vC=[]; vP=[]
for lam in S:
    for mu in S:
        c=solve_c(S,A,lam,mu)
        for nu in S:
            totA+=1
            if (c[nu]!=0) and not (E.geq(nu,lam) and E.geq(nu,mu)): vA.append((lam,mu,nu))
            else: nA+=1
        if E.geq(lam,mu):
            totC+=1
            if sp.simplify(sp.cancel(c[lam]-A[(lam,mu)]))==0: nC+=1
            else: vC.append((lam,mu))
        for nu in S:
            if c[nu]==0: continue
            totP+=1
            if is_laurent(c[nu]): nP+=1
            else: vP.append((lam,mu,nu))
    print(f'  ...lam={lam} done [{time.time()-t0:.0f}s]',flush=True)
print(f'N={N}  Theorem A support bound : {nA}/{totA}   violations {vA[:6]}')
print(f'N={N}  Cor C extremal value    : {nC}/{totC}   violations {vC[:6]}')
print(f'N={N}  (G2) Laurent in w       : {nP}/{totP}')
print(f'N={N}  (G2) failures           : {vP}')
