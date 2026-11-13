from fphase import *
from collections import defaultdict
rel=defaultdict(lambda: defaultdict(set)); n=0
for a in range(1,5):
 for c_ in range(1,5):
  for lam in parts(a,4):
   for mu in parts(c_,4):
    for nu,bs in LR(lam,mu,4).items():
      if len(bs)<2: continue
      nt=trim(nu)
      for bb in bs:
        T=beta_cm(bb,lam,mu,nt); land,seq=ftrace(T,nt); tr=track(T,nt); n+=1
        for x in T: rel[lam][tr[x]].add(land[x])
bad=0
for lam,R in rel.items():
    f=all(len(v)==1 for v in R.values())
    if not f: bad+=1
    print(lam, 'function' if f else 'NOT a function', {k:sorted(v) for k,v in sorted(R.items())})
print('elements',n,'lambdas where endpoint->landing is not a function:',bad,'of',len(rel))
