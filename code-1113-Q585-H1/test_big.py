"""Prop 5.1 (R_ii=R_iii), K1 (kappa = hw.revrot.Tao), kappa(M1)=T, H1 on c>=2 cells |lam|,|mu|<=5 (the 134-element census).
Controls: K0 (revrot without raising) and Klex (hw of revrot, then mapped to the lex-other LR tableau) must fail."""
from test_K import hw, revrot
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
from switch import sigma
import time
st={}; t0=time.time()
def bump(k,c=True): st[k]=st.get(k,0)+bool(c)
for a in range(1,6):
 for c_ in range(1,6):
  for lam in parts(a,4):
   for mu in parts(c_,4):
    for nu,bs in LR(lam,mu,4).items():
      if len(bs)<2: continue
      nt=trim(nu); bump('cells'); cellK1=True
      for bb in bs:
        T=beta_cm(bb,lam,mu,nt)
        k0=max(nt[0],1)+1; k1=len(nt)+1
        b=canonical(nt,k0,k1); B=b.nest[1]; EB,NB=B['u'],B['w']
        F={B['pos'][(c,r)]:v for (c,r),v in T.items()}
        migrate_flock(b,1,F,EB,NB,'E')
        _,Rii=readout(b,0,1,'N'); _,Riii=readout(b,0,-1,'W'); _,K=readout(b,0,-1,'S'); _,PG=readout(b,1,-1,'W')
        bump('n'); bump('P5.1',norm(Rii)==norm(Riii)); bump('kappa=T',norm(K)==norm(T))
        k1ok=norm(hw(revrot(Rii,len(lam))))==norm(T); bump('K1',k1ok); cellK1&=k1ok
        bump('H1',norm(PG)==norm(sigma(T,nt)))
        bump('K0 (control)',norm(revrot(Rii,len(lam)))==norm(T))
      bump('cells K1',cellK1)
print(st,'%.0fs'%(time.time()-t0))
