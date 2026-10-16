"""Q492 FULL VERIFICATION, with pre-registered controls and two plants.

CONTROL 1 (exact pole divisor; independent mechanism = pure combinatorics):
  predicted long-root pole set of c^nu_{lam,mu}
     = { H_ab : exists kappa <= nu, kappa >= lam, kappa >= mu,
                {a,b} <= kappa, kappa\{a,b} >= lam, kappa\{a,b} >= mu }
  (poles are CREATED at kappa and PROPAGATE up the recursion)
  compared against the measured denominator from SymPy.

CONTROL 2 (collision theorem): no measured pole ever lies on a collision
  divisor (w_a = +-1, w_a = w_b) or a q-divisor.

PLANT P1 (K-matrix off-diagonal): (u_j^2 - 1) -> (u_j^2 - q)   [kills the
  short-root collision structure]
PLANT P2 (R-matrix weight): (u_j u_k - 1) -> (u_j u_k - q) in the j<k<=r block
  [kills the long-root factor]
Each plant must be reported FIRED or SPARED; a spared plant is informative.
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp, itertools
from collections import Counter
import q484_engine as E
q=E.q
geq=E.geq

def linext(S): return sorted(S,key=lambda M:(len(M),M))

def predicted_poles(N,S,lam,mu,nu):
    out=set()
    for kap in S:
        if not (geq(nu,kap) and geq(kap,lam) and geq(kap,mu)): continue
        for a in range(1,N+1):
            for b in range(a+1,N+1):
                if a in kap and b in kap:
                    km=tuple(x for x in kap if x not in (a,b))
                    if geq(km,lam) and geq(km,mu): out.add((a,b))
    return out

def measured_poles(f,w,N):
    f=sp.cancel(sp.together(f))
    den=sp.expand(sp.fraction(f)[1])
    lr=set(); other=[]
    for fac,m in sp.factor_list(den)[1]:
        if len(sp.Poly(fac,*w).monoms())==1: continue
        hit=None
        for a in range(1,N+1):
            for b in range(a+1,N+1):
                t=sp.expand(w[a-1]*w[b-1]-1)
                if sp.expand(fac-t)==0 or sp.expand(fac+t)==0: hit=(a,b)
        if hit: lr.add(hit)
        else: other.append(sp.factor(fac))
    return lr,other

def run(N, Fpatch=None, tag='clean'):
    w=E.wvars(N); S=linext(E.subsets(N))
    if Fpatch: E.F=Fpatch
    A={(M,I):E.Atilde(N,w,M,I,'C') for M in S for I in S}
    def solve(lam,mu):
        c={}
        for nu in S:
            v=A[(nu,lam)]*A[(nu,mu)]; acc=sp.Integer(0)
            for kap in S:
                if kap==nu: break
                if c.get(kap,0)!=0: acc+=c[kap]*A[(nu,kap)]
            c[nu]=sp.cancel(sp.together((v-acc)/A[(nu,nu)]))
        return c
    n_div_ok=n_div=0; n_coll_ok=n_coll=0; mism=[]
    for lam in S:
        for mu in S:
            if mu<lam: continue
            c=solve(lam,mu)
            for nu in S:
                if c[nu]==0: continue
                pred=predicted_poles(N,S,lam,mu,nu)
                meas,other=measured_poles(c[nu],w,N)
                n_div+=1
                if pred==meas: n_div_ok+=1
                else: mism.append((lam,mu,nu,sorted(pred),sorted(meas)))
                n_coll+=1
                if not other: n_coll_ok+=1
                else: mism.append((lam,mu,nu,'NONLONGROOT',other))
    print(f'  [{tag}] N={N}  CONTROL 1 exact pole divisor : {n_div_ok}/{n_div}')
    print(f'  [{tag}] N={N}  CONTROL 2 no collision poles : {n_coll_ok}/{n_coll}')
    for m in mism[:6]: print('     mismatch:',m)
    return n_div_ok==n_div and n_coll_ok==n_coll

print('=== CLEAN RUN ===')
import q484_engine as E0
F_orig=E.F
ok2=run(2,tag='clean'); 
print('=== PLANTS at N=2 ===')
from math import factorial
from itertools import permutations
def make_F(plant):
    def F(N,n,r,u,w,I,beta,gamma):
        total=sp.Integer(0)
        for sigma in permutations(range(n)):
            us=[u[sigma[j]] for j in range(n)]; t=sp.Integer(1)
            for j in range(r):
                ij=I[j]
                sq = (us[j]**2 - q) if plant=='P1' else (us[j]**2 - 1)
                t*=(1-q)*sq
                for k in range(ij+1,N+1): t*=(us[j]-q*w[k-1])
                for k in range(1,ij): t*=(us[j]-w[k-1])
            for j in range(r):
                for k in range(j+1,r):
                    pr = (us[j]*us[k]-q) if plant=='P2' else (us[j]*us[k]-1)
                    t*=(q*us[j]-us[k])/(us[j]-us[k])*pr
            for j in range(r,n):
                for k in range(N): t*=(us[j]-q*w[k])
            for j in range(r,n):
                for k in range(r): t*=(q*us[j]-us[k])/(us[j]-us[k])*(us[j]*us[k]-1)
            for j in range(r,n):
                for k in range(j+1,n): t*=(us[j]*us[k]-q)
            for j in range(r,n): t*=(beta*us[j]-gamma)
            total+=t
        pref=sp.Integer(1)
        for k in range(n): pref*=u[k]**(-(N+n))
        return sp.together(pref*total/factorial(n-r))
    return F
for p in ('P1','P2'):
    r=run(2,Fpatch=make_F(p),tag=p)
    print(f'  PLANT {p}: {"SPARED (controls still pass) -- INFORMATIVE" if r else "FIRED (controls broke)"}')
    E.F=F_orig
