"""Uniqueness controls for Lemma D(1), all KT puzzles n<=3 (2, 5, 14; same enumeration as
scratch/q597_k1_emit.py, re-coded here).  For each puzzle:
 (a) the 0/1-edge graph on the vertices of T_n is connected  (=> lift unique up to translation);
 (b) lesion: for each single 0/1-edge e, does deleting it disconnect the graph? (if so a second,
     non-translate lift exists: shift f on one component); count bridges.
 (c) lesion: keep only the 0-edges (drop all 1-edges): connected? (predicted: no, when a 1 exists)
"""
import itertools
def ups(n): return [(a,b) for a in range(n) for b in range(n-a)]
def downs(n): return [(a,b) for a in range(n) for b in range(n-a-1)]
def upE(a,b): return [(a,b,0),(a,b,2),(a,b,1)]
def dnE(a,b): return [(a+1,b,1),(a,b+1,0),(a,b,2)]
KT=[(0,0,0),(1,1,1),(10,0,1),(0,1,10),(1,10,0)]
def boundary(n): return [(i,0,0) for i in range(n)]+[(n-1-j,j,2) for j in range(n)]+[(0,j,1) for j in range(n)]
def ends(e):
    a,b,k=e
    return {0:((a,b),(a+1,b)),1:((a,b),(a,b+1)),2:((a+1,b),(a,b+1))}[k]
def ncomp(n,edges):
    V=[(a,b) for a in range(n+1) for b in range(n+1-a)]
    par={v:v for v in V}
    def f(x):
        while par[x]!=x: x=par[x]
        return x
    for e in edges:
        p,q=ends(e); par[f(p)]=f(q)
    return len({f(v) for v in V})
tot=conn=0; bridges=nb=0; only0_disc=has1=0
for n in (1,2,3):
    for choice in itertools.product(KT, repeat=len(ups(n))):
        lab={}
        for t,p in zip(ups(n),choice): lab.update(zip(upE(*t),p))
        if any(lab[e]==10 for e in boundary(n)): continue
        if not all(tuple(lab[e] for e in dnE(*t)) in KT for t in downs(n)): continue
        tot+=1
        E01=[e for e,x in lab.items() if x!=10]
        conn+= ncomp(n,E01)==1
        for e in E01:
            nb+=1; bridges+= ncomp(n,[g for g in E01 if g!=e])>1
        if any(x==1 for x in lab.values()):
            has1+=1; only0_disc+= ncomp(n,[e for e,x in lab.items() if x==0])>1
print(f"puzzles {tot}; 0/1-graph connected {conn}/{tot}")
print(f"lesion single-edge deletion: {bridges}/{nb} deletions disconnect (=> a non-translate 2nd lift)")
print(f"lesion 0-edges only: disconnected in {only0_disc}/{has1} puzzles containing a 1")
