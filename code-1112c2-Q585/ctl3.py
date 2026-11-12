import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112-Q585')
from commutor import LR, C1, parts
from collections import Counter
NC=4; S=Counter()
for a in range(1,8):
  for c_ in range(1,8):
    for lam in parts(a,NC):
      for mu in parts(c_,NC):
        L1=LR(lam,mu,NC); L2=None
        for nu,bs in L1.items():
          if len(bs)<3: continue
          if L2 is None: L2=LR(mu,lam,NC)
          src=sorted(bs); tgt=sorted(L2[nu]); Ci={b:C1(b,lam,mu,NC)[0] for b in bs}
          S['cells c=%d'%len(bs)]+=1
          S['lex differs']+=any(tgt[src.index(b)]!=Ci[b] for b in bs)
          S['revlex differs']+=any(tgt[::-1][src.index(b)]!=Ci[b] for b in bs)
print(dict(S))
