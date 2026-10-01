"""Verify: beta = c_+ where c(r) = min( H(r)-Lo(r)+1 , N*(Lambda+1-r) ) is CONCAVE.
   H(r)  = min(T+, t0+r, omega-r)      concave
   Lo(r) = max(T-, t0-r, alpha+r)      convex
   N = T+ - T- + 1,  Lambda = sigma0 - max(u, A+B+1),  alpha = tau_wedge - Lambda,
   omega = tau_vee + Lambda."""
from regions import *
from fractions import Fraction as F

def true_beta(n,A,B,u,Tm,Tp):
    beta={}
    for t in range(Tm,Tp+1):
        at,bt,ct,dt=endpoints(t,n,A,B,u)
        nt,mt=bt-at+1,dt-ct+1
        if nt<1 or mt<1: continue
        L=F(nt+mt-2,2); dl=F(abs(nt-mt),2)
        r=dl
        while r<=L: beta[r]=beta.get(r,0)+1; r=r+1
    return beta

def formula_beta(n,A,B,u,Tm,Tp):
    C,D=B+1,A+n-1
    s0=F(u+B+D,2); t0=F(u+B-D,2)
    Lam=s0-max(u,A+B+1)
    tau1,tau2=A,u-1-B
    tw,tv=min(tau1,tau2),max(tau1,tau2)
    alpha=tw-Lam; omega=tv+Lam
    N=Tp-Tm+1
    out={}
    # lattice of r: same class as L_t, i.e. Z - t0  (equivalently r+t0 in Z)
    r = F(0) if t0==int(t0) else F(1,2)
    while r<=Lam+25:
        H=min(F(Tp),t0+r,omega-r); Lo=max(F(Tm),t0-r,alpha+r)
        c=min(H-Lo+1, N*(Lam+1-r))
        if c>0: out[r]=int(c)
        r=r+1
    return out, Lam, alpha, omega, N

def concave_on_lattice(n,A,B,u,Tm,Tp):
    """check c is concave in r on the lattice"""
    C,D=B+1,A+n-1
    s0=F(u+B+D,2); t0=F(u+B-D,2); Lam=s0-max(u,A+B+1)
    tau1,tau2=A,u-1-B; tw,tv=min(tau1,tau2),max(tau1,tau2)
    alpha=tw-Lam; omega=tv+Lam; N=Tp-Tm+1
    cf=lambda r: min(min(F(Tp),t0+r,omega-r)-max(F(Tm),t0-r,alpha+r)+1, N*(Lam+1-r))
    r0 = F(0) if t0==int(t0) else F(1,2)
    r=r0-10
    while r<=Lam+15:
        if 2*cf(r) < cf(r-1)+cf(r+1): return False,r
        r=r+1
    return True,None

bad=0; nonconc=0; tot=0; ex=None; exc=None
for n in range(2,12):
  for mu1 in range(0,2):
    for mu2 in range(mu1,mu1+n+1):
      for lam1 in range(mu1,mu1+13):
        for lam2 in range(max(mu2,lam1),lam1+n+1):
          d=(lam1-mu1)+(lam2-mu2)
          if d>14 or d<1 or not valid_shape(n,(mu1,mu2),(lam1,lam2)): continue
          for b in range(0,d+1):
            nn,A,B,u,Tm,Tp=cyl_params(n,(mu1,mu2),(lam1,lam2),b)
            if Tm>Tp: continue
            tb=true_beta(n,A,B,u,Tm,Tp)
            if not tb: continue
            tot+=1
            fb,Lam,al,om,N=formula_beta(n,A,B,u,Tm,Tp)
            if fb!=tb:
                bad+=1
                if ex is None: ex=(n,(mu1,mu2),(lam1,lam2),b,dict(sorted(tb.items())),dict(sorted(fb.items())))
            ok,r=concave_on_lattice(n,A,B,u,Tm,Tp)
            if not ok:
                nonconc+=1
                if exc is None: exc=(n,(mu1,mu2),(lam1,lam2),b,r)
print(f"instances={tot}")
print(f"formula beta != true beta : {bad}   {ex}")
print(f"c not concave on lattice  : {nonconc}   {exc}")
print()
print("--- control: DROP the N*(Lambda+1-r) clamp; it must now MISMATCH ---")
def formula_noclamp(n,A,B,u,Tm,Tp):
    C,D=B+1,A+n-1
    s0=F(u+B+D,2); t0=F(u+B-D,2); Lam=s0-max(u,A+B+1)
    tau1,tau2=A,u-1-B; tw,tv=min(tau1,tau2),max(tau1,tau2)
    alpha=tw-Lam; omega=tv+Lam
    out={}; r=F(0) if t0==int(t0) else F(1,2)
    while r<=Lam+25:
        c=min(F(Tp),t0+r,omega-r)-max(F(Tm),t0-r,alpha+r)+1
        if c>0: out[r]=int(c)
        r=r+1
    return out
bad2=0; tot2=0; ex2=None
for n in range(2,12):
  for mu2 in range(0,n+1):
    for lam1 in range(0,13):
      for lam2 in range(max(mu2,lam1),lam1+n+1):
        d=lam1+(lam2-mu2)
        if d>14 or d<1 or not valid_shape(n,(0,mu2),(lam1,lam2)): continue
        for b in range(0,d+1):
          nn,A,B,u,Tm,Tp=cyl_params(n,(0,mu2),(lam1,lam2),b)
          if Tm>Tp: continue
          tb=true_beta(n,A,B,u,Tm,Tp)
          if not tb: continue
          tot2+=1
          if formula_noclamp(n,A,B,u,Tm,Tp)!=tb:
              bad2+=1
              if ex2 is None: ex2=(n,(0,mu2),(lam1,lam2),b)
print(f"  without clamp: mismatches {bad2}/{tot2}  first={ex2}")
assert bad2>0, "CONTROL BLIND: clamp appears unnecessary, re-examine"
