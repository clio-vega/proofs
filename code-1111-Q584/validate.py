"""Global control: unit-triangle {0,1,X} model on side-n triangle == KTW tile count_puzzles."""
from star import *
from puzzles import count_puzzles, boundary_labels, region_tris
import itertools
sys.setrecursionlimit(10000)
def tri_count(n,lam,mu,nu,k,upK=0,dnK=0):
    bl=boundary_labels(n,lam,mu,nu,k)
    tris=[(t,a,b) for t,a,b in region_tris(n)]
    # simple DFS row by row with edge constraints
    L=dict(bl); out={}
    def rec(i,ku,kd):
        if i==len(tris):
            out[(ku,kd)]=out.get((ku,kd),0)+1; return
        t=tris[i]; es=tri_edges(*t); S=UP if t[0]=='U' else DN
        cands=list(S)+([(X,X,X)] if (upK if t[0]=='U' else dnK) else [])
        for trip in cands:
            new=[];bad=False
            for e,l in zip(es,trip):
                if e in L:
                    if L[e]!=l: bad=True;break
                else: new.append(e)
            if bad: continue
            for e,l in zip(es,trip):
                if e not in L: L[e]=l
            isK = trip==(X,X,X) and trip not in S
            rec(i+1,ku+(isK and t[0]=='U'),kd+(isK and t[0]=='D'))
            for e in new: del L[e]
    rec(0,0,0); return out
def parts(k,m):
    return [p for p in itertools.product(range(m,-1,-1),repeat=k) if all(p[i]>=p[i+1] for i in range(k-1))]
bad=tot=0
for n,k in ((3,1),(4,2),(5,2)):
    P=parts(k,n-k)
    for la,mu,nu in itertools.product(P,P,P):
        a=count_puzzles(n,list(la),list(mu),list(nu),k,['zero_A','zero_B','one_A','one_B','rh_014','rh_230','rh_125'])
        b=tri_count(n,list(la),list(mu),list(nu),k).get((0,0),0)
        tot+=1; bad+=(a!=b)
        if a!=b: print('MISMATCH',n,la,mu,nu,a,b)
print('unit-triangle vs KTW tiles: %d/%d agree'%(tot-bad,tot))
