"""Q584 V2: cyclic class order of bulk legs for regions with exactly 4 bulk legs.
Question: is every 4-bulk-leg region 'rhombus-like' (opposite bulk legs equal class: x,y,x,y)?"""
from regions import census, tri_dedges
from puzzles import D, edge_key
from collections import Counter, defaultdict
import sys
S=int(sys.argv[1]) if len(sys.argv)>1 else 9
allR,*_=census(S)
pat=Counter(); skipped=0; examples={}
for R in allR:
    cnt=Counter()
    for t in R:
        for p,d in tri_dedges(t): cnt[edge_key(p,d)]+=1
    bd=[(p,d) for t in R for p,d in tri_dedges(t) if cnt[edge_key(p,d)]==1]
    start=defaultdict(list)
    for p,d in bd: start[p].append(d)
    if any(len(v)>1 for v in start.values()): skipped+=1; continue   # pinch point
    # walk
    p0,d0=bd[0]; seq=[]; p,d=p0,d0; comps=1; used=0
    while True:
        seq.append((p,d)); used+=1
        q=(p[0]+D[d][0],p[1]+D[d][1]); d=start[q][0]; p=q
        if (p,d)==(p0,d0): break
    if used!=len(bd): skipped+=1; continue   # multiple boundary components (holes)
    bulk=[d%3 for p,d in seq if not (d%3==1 and p[0]==0)]
    if len(bulk)!=4: continue
    # canonical cyclic pattern
    rots=[tuple(bulk[i:]+bulk[:i]) for i in range(4)]
    c=min(rots); kind='opposite-equal' if c[0]==c[2] and c[1]==c[3] else 'ADJACENT'
    pat[(c,kind,any(d%3==1 and p[0]==0 for p,d in seq))]+=1
    examples.setdefault((c,kind),sorted(R))
print('skipped (pinch/holes):',skipped)
for k,v in sorted(pat.items()): print(k,v)
for k,v in examples.items():
    if k[1]=='ADJACENT': print('example',k,v)
