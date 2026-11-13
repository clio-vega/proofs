"""Per G journey: ordered hexagon-step types vs the moves of the paired Sigma slide (slide |mu|+1-k)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from trace1 import *
import migrate

def sigma_moves(S, outer):
    """switch.sigma, recording each slide's move string ('U' next row, 'R' right) and its vacated cell"""
    S = dict(S); outer=list(outer); inner=[0]*len(outer)
    for r,l in enumerate(outer):
        cs=[c for (c,rr) in S if rr==r]; inner[r]=min(cs) if cs else l
    kappa=[x for x in inner if x]
    cells=[(c,r) for r,l in enumerate(kappa) for c in range(l)]
    cells.sort(key=lambda x:(-x[1],-x[0]))
    res=[]
    for (c,r) in cells:
        hc,hr=c,r; mv=''
        while True:
            R=S.get((hc+1,hr)); D=S.get((hc,hr+1))
            if R is None and D is None: break
            if D is not None and (R is None or D<=R):
                S[(hc,hr)]=D; del S[(hc,hr+1)]; hr+=1; mv+='U'
            else:
                S[(hc,hr)]=R; del S[(hc+1,hr)]; hc+=1; mv+='R'
        res.append(((c,r),(hc,hr),r+1,mv))
    return res

# monkeypatch journey to record ordered step types
orig_step = migrate.Board.step
LOG=[]
def step(self,R,d):
    out=orig_step(self,R,d); LOG.append(self.last); return out
migrate.Board.step=step

def gtrace(T,nt):
    global LOG
    k0=max(nt[0],1)+1; k1=len(nt)+1
    b=canonical(nt,k0,k1); B=b.nest[1]; EB,NB=B['u'],B['w']; neg=lambda v: tuple(-a for a in v)
    F={B['pos'][(c,r)]:v for (c,r),v in T.items()}
    mu_pos=[ij for ij in B['pos'] if ij not in T]
    boxesG={box_of(b,1,b.tiles[B['pos'][ij]][1],neg(EB),neg(NB)):ij for ij in mu_pos}
    mu_shape=tuple(sorted(Counter(r for _,r in mu_pos).values(),reverse=True))
    Ls=lr_tableaux(set(boxesG),mu_shape)
    Gfl={B['pos'][boxesG[bx]]:v for bx,v in Ls[0].items()}
    traced(b,1,F,EB,NB,'E')
    # G phase, journey by journey
    d=unit_dir(EB)  # W of (-E) = E_B direction
    E_,N_=neg(EB),neg(NB)
    key={t:(v,box_of(b,1,b.tiles[t][1],E_,N_)[0]) for t,v in Gfl.items()}
    order=sorted(Gfl,key=lambda t:key[t])
    inv={t:ij for ij,t in b.nest[1]['pos'].items()}
    out=[]
    for t in order:
        del b.nest[1]['pos'][inv[t]]; LOG=[]
        R,X,ij,st=b.journey(t,1,d)
        out.append((Gfl[t],''.join(LOG),R))
    E2,N2,_=induced(b,X,E_,N_)
    return [(v,s,box_of(b,X,b.tiles[R][1],E2,N2)) for v,s,R in out]

if __name__=='__main__':
    import itertools
    n=0
    for a in range(1,5):
     for c_ in range(1,5):
      for lam in parts(a,4):
       for mu in parts(c_,4):
        for nu,bs in LR(lam,mu,4).items():
          if len(bs)<2: continue
          nt=trim(nu)
          for bb in bs:
            T=beta_cm(bb,lam,mu,nt); g=gtrace(T,nt); sm=sigma_moves(T,nt)[::-1]
            n+=1
            if n<=6:
              print('lam',lam,'mu',mu,'nu',nt)
              for (v,s,cell),(start,end,lab,mv) in zip(g,sm):
                  print('   G val',v,'lands',cell,'steps',s,'| slide into',start,'->',end,'lab',lab,'moves',mv)
