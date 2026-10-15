"""Q484 engine plants: perturb ONE factor of eq (Motherfunction) and check that
the diagonal root-product control (Control 1) fails.  A control that cannot be
made to fail is decoration."""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp
from itertools import permutations
from math import factorial
import q484_engine as E
q = E.q

def F_planted(N,n,r,u,w,I,beta,gamma,plant):
    total = sp.Integer(0)
    for sigma in permutations(range(n)):
        us=[u[sigma[j]] for j in range(n)]; t=sp.Integer(1)
        for j in range(r):
            ij=I[j]
            t *= (1-q)*(us[j]**2 - (0 if plant=='P2_drop_u2m1' else 1))
            for k in range(ij+1,N+1): t *= (us[j]-q*w[k-1])
            for k in range(1,ij):     t *= (us[j]-w[k-1])
        for j in range(r):
            for k in range(j+1,r):
                t *= (q*us[j]-us[k])/(us[j]-us[k]) * (us[j]*us[k]-1)
        for j in range(r,n):
            for k in range(N): t *= (us[j]-q*w[k])
        for j in range(r,n):
            for k in range(r):
                t *= (q*us[j]-us[k])/(us[j]-us[k]) * (us[j]*us[k]-1)
        for j in range(r,n):
            for k in range(j+1,n):
                # PLANT P1: the reflecting-boundary pair factor (u_j u_k - q) -> (u_j u_k - 1)
                t *= (us[j]*us[k]-(1 if plant=='P1_pairfactor' else q))
        for j in range(r,n):
            t *= (beta*us[j]-gamma)
        total += t
    pref=sp.Integer(1)
    for k in range(n): pref *= u[k]**(-(N+n))
    return sp.together(pref*total/factorial(n-r))

def Atilde_planted(N,w,M,I,plant):
    s,r=len(M),len(I)
    if s<r: return sp.Integer(0)
    n=s
    if n==0: return sp.cancel(sp.Integer(1)/E.EC(N,0,0,[],w))
    u=[w[m-1] for m in M]
    return sp.cancel(sp.together(F_planted(N,n,r,u,w,I,sp.Integer(1),sp.Integer(0),plant)/E.EC(N,n,r,u,w)))

if __name__=='__main__':
    for plant in ('P1_pairfactor','P2_drop_u2m1'):
        for N in (2,3):
            w=E.wvars(N); bad=0; tot=0
            for I in E.subsets(N):
                a=Atilde_planted(N,w,I,I,plant); d=E.diag_rootproduct(N,w,I)
                tot+=1
                if sp.simplify(sp.cancel(a-d))!=0: bad+=1
            print(f'PLANT {plant}  N={N}: diagonal control fails on {bad}/{tot} labels'
                  f'  -> {"FIRES" if bad else "DID NOT FIRE"}')
