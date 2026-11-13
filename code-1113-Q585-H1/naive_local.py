"""Refutation census for the naive local lemma: is a G journey's hexagon-step string a function of
its paired Sigma slide's move string?  Also: does the pairing journey k <-> slide |mu|+1-k hold (landing cell)?"""
from trace2 import *
from collections import defaultdict
pairs=defaultdict(set); n=land=0; lens=[]
for a in range(1,5):
 for c_ in range(1,5):
  for lam in parts(a,4):
   for mu in parts(c_,4):
    for nu,bs in LR(lam,mu,4).items():
      if len(bs)<2: continue
      nt=trim(nu)
      for bb in bs:
        T=beta_cm(bb,lam,mu,nt); g=gtrace(T,nt); sm=sigma_moves(T,nt)[::-1]
        for (v,s,cell),(start,end,lab,mv) in zip(g,sm):
            n+=1; land+= (cell==end and v==lab); pairs[mv].add(s); lens.append((len(mv),s.count('a'),s.count('b')))
print('journeys',n,'landing matches paired slide (cell & label):',land)
print('distinct move strings',len(pairs),'; step strings per move string:',{k:len(v) for k,v in sorted(pairs.items())})
inj=sum(1 for v in pairs.values() if len(v)==1); print('move strings with a unique step string:',inj)
# is any step string shared by two different move strings?
inv=defaultdict(set)
for k,v in pairs.items():
    for s in v: inv[s].add(k)
print('step strings shared across different move strings:',sum(1 for v in inv.values() if len(v)>1),'of',len(inv))
