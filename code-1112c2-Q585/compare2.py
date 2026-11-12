"""M (Purbhoo Cor 3.4) vs C (HK commutor) under TWO bijections, on all cells with c>=2 in a range; reports by c.
Also: Prop 5.4 involution, Thm 1 (LR output), lex control (must differ from C to be live)."""
import sys, time
from collections import Counter, defaultdict
from compare import M, beta_cm, box, trim, LR, C1, C2, parts, NC as _
from betaq import beta_q
from flock import is_lr, content
NC=int(sys.argv[1]); amax=int(sys.argv[2]); cmin=int(sys.argv[3])
S=defaultdict(Counter); t0=time.time(); witnesses=[]
for a in range(1,amax+1):
  for c_ in range(1,amax+1):
    for lam in parts(a,NC):
      for mu in parts(c_,NC):
        L1=LR(lam,mu,NC); L2=None
        for nu,bs in L1.items():
          c=len(bs)
          if c<cmin: continue
          if L2 is None: L2=LR(mu,lam,NC)
          nu_t=trim(nu); s=S[c]; s['cells']+=1
          Cimg={b:C1(b,lam,mu,NC)[0] for b in bs}
          assert all(Cimg[b]==C2(b,lam,mu,NC)[0] for b in bs)
          cellok={'Q':True,'CM':True}
          for b in bs:
            s['elems']+=1
            for name,fwd,bwd,inner_in,inner_out in (
                ('Q', lambda x: beta_q(x,lam,mu), lambda y: beta_q(y,mu,lam), lam, mu),
                ('CM',lambda x: beta_cm(x,lam,mu,nu_t), lambda y: beta_cm(y,mu,lam,nu_t), mu, lam)):
              T=fwd(b); out=M(T,nu_t)
              if not is_lr(out) or set(out)!=box(nu_t,inner_out): s['bad_output_'+name]+=1
              if M(out,nu_t)!=T: s['invol_bad_'+name]+=1
              ok= out==bwd(Cimg[b]); s['agree_'+name]+=ok; cellok[name]&=ok
              if not ok and len(witnesses)<5: witnesses.append((name,lam,mu,nu_t,b,Cimg[b],T,out,bwd(Cimg[b])))
          for k in cellok: s['cells_agree_'+k]+=cellok[k]
          src=sorted(bs); tgt=sorted(L2[nu])
          s['lex_differs_C']+= any(tgt[src.index(b)]!=Cimg[b] for b in bs)
          s['C_identity_like']+= all(Cimg[b]==b for b in bs)
print('gl%d |lam|,|mu|<=%d, c>=%d  (%.0fs)'%(NC,amax,cmin,time.time()-t0))
for c in sorted(S): print(' c=%d'%c,dict(S[c]))
for w in witnesses: print('WITNESS',w)
