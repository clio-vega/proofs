from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
from switch import sigma
def rot180(T):
    return {(-c,-r):v for (c,r),v in T.items()}
def transp(T): return {(r,c):v for (c,r),v in T.items()}
st={}; ex=[]
def bump(k): st[k]=st.get(k,0)+1
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
        migrate_flock(b,1,F,EB,NB,'E')
        _,Rii=readout(b,0,1,'N'); _,Riii=readout(b,0,-1,'W'); _,K=readout(b,0,-1,'S')
        _,PhiG=readout(b,1,-1,'W'); S=sigma(T,nt)
        bump('n')
        bump('H1') if norm(PhiG)==norm(S) else None
        for nm,f in (('id',lambda x:x),('rot',rot180),('tr',transp),('rot.tr',lambda x:rot180(transp(x)))):
            if norm(f(Rii))==norm(Riii): bump('Rii=Riii via '+nm)
            if norm(f(Rii))==norm(K): bump('Rii=kappa via '+nm)
            if norm(f(Rii))==norm(S): bump('Rii=Sigma(T) via '+nm)
            if norm(f(Riii))==norm(S): bump('Riii=Sigma(T) via '+nm)
        if len(ex)<3: ex.append((lam,mu,nt,norm(T),norm(Rii),norm(Riii),norm(S)))
print(st)
for e in ex: print(e)
