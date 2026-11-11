"""Q584 V1: exhaustive polyiamond census in the half-plane a>=0 (wall a=0).
For every union R of unit triangles (connected, size<=S, anchored to touch a bounded window),
boundary edges split into WALL edges (on a=0, class 1) and BULK edges.  Record bulk class multiset.
Checks:  (L1) every region with exactly 2 bulk edges has bulk classes {0,2} and touches the wall;
         (L2) no region has all bulk edges in one class c (wall edges allowed);
         (L3) closure n0=n2=-n1 for every region (signed counts, ccw orientation) -- instrument check.
Positive controls: the wall triangle U(0,b) must appear with bulk {0,2}; a rhombus with {0,0,1,1}."""
import sys, itertools
from collections import Counter
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import D, edge_key
def tri_dedges(t):
    kind,a,b=t
    P=[(a,b),(a+1,b),(a,b+1)] if kind=='U' else [(a+1,b),(a+1,b+1),(a,b+1)]
    out=[]
    for i in range(3):
        p,q=P[i],P[(i+1)%3]; d=D.index((q[0]-p[0],q[1]-p[1])); out.append((p,d))
    return out
def nbrs(t):
    kind,a,b=t
    if kind=='U': return [('D',a,b-1),('D',a-1,b),('D',a,b)]
    return [('U',a,b),('U',a+1,b),('U',a,b+1)]
def valid(t):  # inside half-plane a>=0: all vertices a>=0
    kind,a,b=t
    return a>=0 if kind=='U' else a>=0   # D(a,b) has vertices a, a+1 -> a>=0
def census(S, B=6):
    seen=set(); stats=Counter(); bad=[]; two=Counter(); mono=0; closure_bad=0
    start=[('U',0,0),('D',0,0),('U',1,0),('D',1,0)]
    frontier=[frozenset([t]) for t in start]
    for s in frontier: seen.add(s)
    level=list(frontier)
    allR=list(frontier)
    for size in range(2,S+1):
        nxt=[]
        for R in level:
            for t in R:
                for u in nbrs(t):
                    if valid(u) and u not in R and -B<=u[2]<=B and u[1]<=B:
                        R2=R|{u}
                        # normalise b-translation (wall is b-invariant): shift min b to 0
                        mb=min(x[2] for x in R2); R2=frozenset((k,a,b-mb) for k,a,b in R2)
                        if R2 not in seen: seen.add(R2); nxt.append(R2)
        level=nxt; allR+=nxt
    for R in allR:
        cnt=Counter()
        for t in R:
            for p,d in tri_dedges(t): cnt[edge_key(p,d)]+=1
        bd=[]
        for t in R:
            for p,d in tri_dedges(t):
                if cnt[edge_key(p,d)]==1: bd.append((p,d))
        n=Counter()
        for p,d in bd: n[d%3]+= (1 if d<3 else -1)
        if not (n[0]==n[2]==-n[1]): closure_bad+=1
        wall=[(p,d) for p,d in bd if d%3==1 and p[0]==0]
        bulk=[(p,d) for p,d in bd if not (d%3==1 and p[0]==0)]
        cls=tuple(sorted(d%3 for p,d in bulk))
        if len(bulk)==2:
            two[(cls,len(wall)>0)]+=1
        if len(set(cls))==1: mono+=1; bad.append((sorted(R),cls))
        stats[len(R)]+=1
    return allR,stats,two,mono,bad,closure_bad
if __name__=='__main__':
    S=int(sys.argv[1]) if len(sys.argv)>1 else 8
    allR,stats,two,mono,bad,cb=census(S)
    print('regions by size:',dict(stats),' total',len(allR))
    print('(L3) closure violations:',cb)
    print('(L1) regions with exactly 2 bulk edges, by (classes, touches wall):',dict(two))
    print('(L2) regions with single-class bulk boundary:',mono, bad[:3])
    # positive controls
    wt=frozenset([('U',0,0)]); print('control wall triangle present:',wt in set(allR))
    trap=frozenset([('U',0,0),('D',0,0),('U',0,1)])
    print('trap region present:',trap in set(allR))
