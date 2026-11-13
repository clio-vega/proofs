"""(a) kappa(M1) == T  [Thm 2 positive control];  also which (sgn,dir) on A reach B / C."""
from readouts import *
from trace1 import run
from compare import beta_cm, trim
from commutor import LR, parts
from switch import sigma
import itertools
n=ok=0; dirs={}
for a in range(1,5):
 for c_ in range(1,5):
  for lam in parts(a,4):
   for mu in parts(c_,4):
    for nu,bs in LR(lam,mu,4).items():
      if len(bs)<2: continue
      nt=trim(nu)
      for bb in bs:
        T=beta_cm(bb,lam,mu,nt)
        k0=max(nt[0],1)+1; k1=len(nt)+1
        b=canonical(nt,k0,k1); B=b.nest[1]; EB,NB=B['u'],B['w']
        F={B['pos'][(c,r)]:v for (c,r),v in T.items()}
        migrate_flock(b,1,F,EB,NB,'E')          # b is now M1 (nest B still holds mu, unlabelled)
        # kappa needs nest B to be just beta; it is (mu). read A:
        for sgn in (1,-1):
          for d in 'NESW':
            try:
              dst,tab=readout(b,0,sgn,d)
              dirs.setdefault((sgn,d),set()).add(dst)
            except Exception as e: dirs.setdefault((sgn,d),set()).add('ERR')
        dst,tab=readout(b,0,-1,'S'); n+=1; ok+= (dst==1 and norm(tab)==norm(T))
print('kappa(M1)==T:',ok,'/',n); print('reach:',dirs)
