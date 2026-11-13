"""Do kappa (-E,-N; south) and R_ii (E,N; north) move each A-rhombus to the same B cell, and end in the same mosaic?"""
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
def journeys(board,sgn,dirn):
    b=copy.deepcopy(board); u,w=b.nest[0]['u'],b.nest[0]['w']
    E_=tuple(sgn*a for a in u); N_=tuple(sgn*a for a in w); fl=unique_flock(b,0,E_,N_)
    src={t:ij for ij,t in b.nest[0]['pos'].items()}
    dvec={'E':E_,'N':N_,'W':tuple(-a for a in E_),'S':tuple(-a for a in N_)}[dirn]; d=unit_dir(dvec)
    key={t:(v,box_of(b,0,b.tiles[t][1],E_,N_)[0]) for t,v in fl.items()}
    order=sorted(fl,key=lambda t:key[t],reverse=(dirn in 'NE')); res={}
    for t in order:
        del b.nest[0]['pos'][src[t]]; R,X,ij,st=b.journey(t,0,d); res[src[t]]=ij
    tiles=frozenset((k,frozenset(vs)) for k,vs in b.tiles.values())
    return res,tiles,[src[t] for t in order]
st={}
def bump(k,c=True): st[k]=st.get(k,0)+bool(c)
for a in range(1,5):
 for c_ in range(1,5):
  for lam in parts(a,4):
   for mu in parts(c_,4):
    for nu,bs in LR(lam,mu,4).items():
      if len(bs)<2: continue
      nt=trim(nu)
      for bb in bs:
        T=beta_cm(bb,lam,mu,nt); k0=max(nt[0],1)+1; k1=len(nt)+1
        b=canonical(nt,k0,k1); B=b.nest[1]; F={B['pos'][(c,r)]:v for (c,r),v in T.items()}
        migrate_flock(b,1,F,B['u'],B['w'],'E')
        j1,t1,o1=journeys(b,-1,'S'); j2,t2,o2=journeys(b,1,'N')
        bump('n'); bump('same final mosaic',t1==t2); bump('same cellwise journeys',j1==j2); bump('same order',o1==o2); bump('reverse order',o1==o2[::-1])
print(st)
