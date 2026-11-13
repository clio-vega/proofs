"""Out-of-sample: H1 (M == Sigma on beta_CM inputs) on gl_4 cells with c>=2 and max(|lam|,|mu|)==6 (new vs c2 census),
plus controls: tie->right lesion of Sigma and lex map must differ."""
import sys, time; sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from switch import sigma
from commutor import LR, C1, parts
from flock import is_lr, content
from compare import beta_cm, M, trim
import migrate
t0=time.time(); st=dict(cells=0,elems=0,H1=0,H2CM=0,lesion_cells=0,lex_cells=0,ties_met=0,c_hist={})
for a in range(1,7):
  for bb in range(1,7):
    if max(a,bb)<6: continue
    for lam in parts(a,4):
      for mu in parts(bb,4):
        for nu,bs in LR(lam,mu,4).items():
          if len(bs)<2: continue
          nt=trim(nu); st['cells']+=1; st['c_hist'][len(bs)]=st['c_hist'].get(len(bs),0)+1
          L2=sorted(LR(mu,lam,4)[nu]); src=sorted(bs); les=lex=False
          for b in bs:
            st['elems']+=1; T=beta_cm(b,lam,mu,nt); s=sigma(T,nt)
            m=M(T,nt); st['H1']+= (m==s)
            cb=C1(b,lam,mu,4)[0]; st['H2CM']+= (s==beta_cm(cb,mu,lam,nt))
            les|= sigma(T,nt,tie='right')!=s; lex|= beta_cm(L2[src.index(b)],mu,lam,nt)!=s
          st['lesion_cells']+=les; st['lex_cells']+=lex
          print(lam,mu,nt,len(bs),'H1 so far',st['H1'],'/',st['elems'],'%.0fs'%(time.time()-t0),flush=True)
st['ties_met']=migrate.Global_ties[:]
print('FINAL',st,'time %.0fs'%(time.time()-t0))
