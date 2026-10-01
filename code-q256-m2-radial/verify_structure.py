from regions import *
from structure import *

print("=== A. structural identities on a wide cylindric sweep ===")
bad=0; checked=0; firstbad=None
for n in range(2,13):
  for mu1 in range(0,n+1):
    for mu2 in range(mu1,mu1+n+1):
      for lam1 in range(mu1,mu1+14):
        for lam2 in range(max(mu2,lam1),lam1+n+1):
          d=(lam1-mu1)+(lam2-mu2)
          if d>14 or d<1 or not valid_shape(n,(mu1,mu2),(lam1,lam2)): continue
          for b in range(0,d+1):
            nn,A,B,u,Tm,Tp = cyl_params(n,(mu1,mu2),(lam1,lam2),b)
            if Tm>Tp: continue
            checked+=1
            res=check_identities(n,A,B,u,Tm,Tp)
            if res:
                bad+=1
                if firstbad is None: firstbad=(n,(mu1,mu2),(lam1,lam2),b,res)
print(f"  slices checked={checked}  identity violations={bad}  {firstbad}")

print()
print("=== B. radial/tail decomposition reproduces G, beta truncated-concave, G log-concave ===")
nG=0; mismatch=0; betabad=0; lcbad=0; firstmm=None; firstbb=None
for n in range(2,13):
  for mu1 in range(0,n+1):
    for mu2 in range(mu1,mu1+n+1):
      for lam1 in range(mu1,mu1+14):
        for lam2 in range(max(mu2,lam1),lam1+n+1):
          d=(lam1-mu1)+(lam2-mu2)
          if d>14 or d<1 or not valid_shape(n,(mu1,mu2),(lam1,lam2)): continue
          for b in range(0,d+1):
            nn,A,B,u,Tm,Tp=cyl_params(n,(mu1,mu2),(lam1,lam2),b)
            if Tm>Tp: continue
            # direct G
            tot={}
            for t in range(Tm,Tp+1):
                for s,v in conv_interval(*endpoints(t,n,A,B,u)).items(): tot[s]=tot.get(s,0)+v
            if not tot: continue
            nG+=1
            s0,beta,gamma,Gr = radial(n,A,B,u,Tm,Tp)
            if Gr!={k:v for k,v in tot.items() if v>0}:
                mismatch+=1
                if firstmm is None: firstmm=(n,(mu1,mu2),(lam1,lam2),b,as_seq(tot),Gr)
            if not beta_is_truncated_concave(beta):
                betabad+=1
                if firstbb is None: firstbb=(n,(mu1,mu2),(lam1,lam2),b,beta)
            if not is_pf2(tot): lcbad+=1
print(f"  nonzero G's={nG}  G!=radial: {mismatch} {firstmm}")
print(f"  beta not truncated-concave: {betabad} {firstbb}")
print(f"  G not PF2 (condition (A) failures): {lcbad}")
