"""F phase: does the F journey of the T-letter at cell x land (nest-A ideal position) where that letter ends up in the
switching Sigma (letters tracked through every jdt move)?  Test both identifications (ij) and (ij transposed)."""
from trace2 import *
import migrate

def track(S, outer):
    """Sigma with letter identities: returns dict start_cell -> final cell (in kappa-region, i.e. u_lambda's cells)"""
    S={c:(v,c) for c,v in S.items()}; outer=list(outer); inner=[0]*len(outer)
    for r,l in enumerate(outer):
        cs=[c for (c,rr) in S if rr==r]; inner[r]=min(cs) if cs else l
    kappa=[x for x in inner if x]; cells=[(c,r) for r,l in enumerate(kappa) for c in range(l)]
    cells.sort(key=lambda x:(-x[1],-x[0]))
    for (c,r) in cells:
        hc,hr=c,r
        while True:
            R=S.get((hc+1,hr)); D=S.get((hc,hr+1))
            if R is None and D is None: break
            if D is not None and (R is None or D[0]<=R[0]): S[(hc,hr)]=D; del S[(hc,hr+1)]; hr+=1
            else: S[(hc,hr)]=R; del S[(hc+1,hr)]; hc+=1
    return {orig:cell for cell,(v,orig) in S.items()}

def ftrace(T,nt):
    k0=max(nt[0],1)+1; k1=len(nt)+1
    b=canonical(nt,k0,k1); B=b.nest[1]; EB,NB=B['u'],B['w']
    F={B['pos'][(c,r)]:(c,r) for (c,r) in T}
    d=unit_dir(EB)
    key={t:(T[F[t]],box_of(b,1,b.tiles[t][1],EB,NB)[0]) for t in F}
    order=sorted(F,key=lambda t:key[t],reverse=True)
    inv={t:ij for ij,t in b.nest[1]['pos'].items()}; out={}; seq=[]
    for t in order:
        del b.nest[1]['pos'][inv[t]]; R,X,ij,st=b.journey(t,1,d); out[F[t]]=ij; seq.append(F[t])
    return out,seq

if __name__=='__main__':
    n=okA=okB=0; ex=None
    for a in range(1,5):
     for c_ in range(1,5):
      for lam in parts(a,4):
       for mu in parts(c_,4):
        for nu,bs in LR(lam,mu,4).items():
          if len(bs)<2: continue
          nt=trim(nu)
          for bb in bs:
            T=beta_cm(bb,lam,mu,nt); land,seq=ftrace(T,nt); tr=track(T,nt)
            n+=1
            okA+= all(land[x]==tr[x] for x in T)
            okB+= all(land[x]==tr[x][::-1] for x in T)
            if ex is None: ex=(lam,mu,nt,T,land,tr,seq)
    print('elements',n,'F landing == switching letter endpoint (ij):',okA,' (transposed):',okB)
    print(ex)
