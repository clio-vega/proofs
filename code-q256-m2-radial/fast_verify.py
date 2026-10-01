from regions import *
from fractions import Fraction as F
from structure import data, check_identities, radial, beta_is_truncated_concave

def beta_pf2(beta):
    if not beta: return True
    ks=sorted(beta); seq=[]; r=ks[0]
    while r<=ks[-1]: seq.append(beta.get(r,0)); r=r+1
    nz=[i for i,v in enumerate(seq) if v>0]
    if nz and nz[-1]-nz[0]+1!=len(nz): return False
    return all(seq[i]**2>=seq[i-1]*seq[i+1] for i in range(1,len(seq)-1))

def Lambda_binding(n,A,B,u,Tm,Tp):
    """Is the r<=Lambda clause ever ACTIVE, i.e. does the 3-interval formula
    give a positive count at some r > Lambda?"""
    C,D,s0,t0,rows=data(n,A,B,u,Tm,Tp)
    Lf=lambda t: F((min(u-1-t,B)-max(t,A)+1)+(min(t+n-1,D)-max(u-t,C)+1)-2,2)
    tau1,tau2=A,u-1-B
    Lam=Lf(min(tau1,tau2))
    assert Lf(max(tau1,tau2))==Lam, "flat piece not flat"
    # alpha, omega from the slope +1 / -1 pieces
    alpha = min(tau1,tau2)-Lam
    omega = max(tau1,tau2)+Lam
    r=Lam+1; hits=0
    while r<=Lam+20:
        hi=min(F(Tp),t0+r,omega-r); lo=max(F(Tm),t0-r,alpha+r)
        if hi-lo+1>0: hits+=1
        r=r+1
    return hits, Lam

tot=0; idbad=0; mm=0; bbad=0; lcbad=0; lambda_active=0
ex_id=ex_mm=ex_b=None
for n in range(2,12):
  for mu1 in range(0,2):                      # translation-invariant in mu1: fix 0,1
    for mu2 in range(mu1,mu1+n+1):
      for lam1 in range(mu1,mu1+13):
        for lam2 in range(max(mu2,lam1),lam1+n+1):
          d=(lam1-mu1)+(lam2-mu2)
          if d>14 or d<1 or not valid_shape(n,(mu1,mu2),(lam1,lam2)): continue
          for b in range(0,d+1):
            nn,A,B,u,Tm,Tp=cyl_params(n,(mu1,mu2),(lam1,lam2),b)
            if Tm>Tp: continue
            G={}
            for t in range(Tm,Tp+1):
                for s,v in conv_interval(*endpoints(t,n,A,B,u)).items(): G[s]=G.get(s,0)+v
            G={k:v for k,v in G.items() if v>0}
            if not G: continue
            tot+=1
            r=check_identities(n,A,B,u,Tm,Tp)
            if r: idbad+=1; ex_id=ex_id or (n,(mu1,mu2),(lam1,lam2),b,r)
            s0,beta,gamma,Gr=radial(n,A,B,u,Tm,Tp)
            if Gr!=G: mm+=1; ex_mm=ex_mm or (n,(mu1,mu2),(lam1,lam2),b,as_seq(G),Gr)
            if not beta_pf2(beta): bbad+=1; ex_b=ex_b or (n,(mu1,mu2),(lam1,lam2),b,beta)
            if not is_pf2(G): lcbad+=1
            h,Lam=Lambda_binding(n,A,B,u,Tm,Tp)
            if h: lambda_active+=1
print(f"slices with G!=0      : {tot}")
print(f"identity violations   : {idbad}  {ex_id}")
print(f"G != gamma(|s-sigma0|): {mm}  {ex_mm}")
print(f"beta not PF2          : {bbad}  {ex_b}")
print(f"G not log-concave  (A): {lcbad}")
print(f"slices where the r<=Lambda clause is ACTIVE (ghost rows): {lambda_active}")
