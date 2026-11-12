"""Q584 Thm 8 under Rick's S1 fix: arcs are chains of rhombi (strand passes through
opposite edges), not single shared edges.  Two wall K's U(0,0) (B) and U(0,bA) (A).
Pattern searched = the RE side signature: one half of A crosses each half of B exactly once
(crossing = a rhombus shared by both strands, using complementary edge pairs), the other half of
A crosses nothing; all other triangles pairwise disjoint; open ends free (triangle beyond unused).
mode 'strict' : strands go straight through every rhombus (Def 7 + S1 fix).
mode 'lesion' : a strand may TURN (adjacent edges) in at most one NON-crossing rhombus.
Prediction: strict 0 (pigeonhole), lesion > 0 (so the 0 is the class argument, not adjacency)."""
import sys, itertools
from collections import Counter
D=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def add(p,q): return (p[0]+q[0],p[1]+q[1])
def U(a,b): return frozenset([(a,b),(a+1,b),(a,b+1)])
def Dn(a,b): return frozenset([(a+1,b),(a,b+1),(a+1,b+1)])
def inH(T): return all(p[0]>=0 for p in T)
def edges(T):
    T=list(T); return [frozenset(e) for e in itertools.combinations(T,2)]
def cls(e):
    p,q=tuple(e); d=(q[0]-p[0],q[1]-p[1]); return D.index(d)%3 if d in D else D.index((-d[0],-d[1]))%3
def other(T,e):
    """the other unit triangle containing edge e"""
    x=(T-e); (x,)=x; p,q=tuple(e); y=(p[0]+q[0]-x[0],p[1]+q[1]-x[1]); return frozenset([p,q,y])
def steps(e,T,turn):
    """strand enters triangle T through edge e.  yield (rhombus=(T,T2), exit edge, turned)."""
    for f in edges(T):
        if f==e: continue
        T2=other(T,f)
        for g in edges(T2):
            if g==f: continue
            straight=(cls(g)==cls(e))
            if straight or turn: yield (T,T2),g,(not straight)
def paths(leg,Kt,L,lesion):
    """all strand halves starting at K's leg edge, <=L rhombi, <=1 turn if lesion"""
    out=[]
    def rec(e,Tprev,rh,used,nturn):
        out.append((tuple(rh),e,nturn))
        if len(rh)==L: return
        T=other(Tprev,e)
        if not inH(T) or T in used: return
        for (A,B),g,t in steps(e,T,lesion):
            if not inH(B) or B in used or (nturn+t)>1: continue
            rec(g,B,rh+[(A,B,e,g)],used|{A,B},nturn+t)
    rec(leg,Kt,[],{Kt},0); return out
def tris(p): return {t for r in p[0] for t in r[:2]}
def beyond(p,Kt):
    e=p[1]; Tl=p[0][-1][1] if p[0] else Kt; return other(Tl,e)
def run(mode,L,bmax):
    lesion=(mode=='lesion'); hits=0; turnwhere=Counter(); reach=0
    KB=U(0,0); legsB=[e for e in edges(KB) if cls(e)!=1]
    PB=[paths(e,KB,L,lesion) for e in legsB]
    for bA in [b for b in range(-bmax,bmax+1) if b!=0]:
        KA=U(0,bA)
        if KA & KB and len(KA&KB)>2: continue
        legsA=[e for e in edges(KA) if cls(e)!=1]
        PA=[paths(e,KA,L,lesion) for e in legsA]
        for p1 in PB[0]:
          t1=tris(p1)
          if KA in t1 or KB in t1: continue
          for p2 in PB[1]:
            t2=tris(p2)
            if t1&t2 or KA in t2: continue
            for i in (0,1):
              for a1 in PA[i]:
                ta=tris(a1)
                if KB in ta or KA in ta: continue
                rh_a={(r[0],r[1]):r for r in a1[0]}
                def cross(pb):
                    c=0
                    for r in pb[0]:
                        k=(r[0],r[1]); k2=(r[1],r[0])
                        ra=rh_a.get(k) or rh_a.get(k2)
                        if ra is None:
                            if (r[0] in ta) or (r[1] in ta): return -1
                            continue
                        # crossing: both straight, complementary edge pairs
                        if cls(ra[2])!=cls(ra[3]) or cls(r[2])!=cls(r[3]): return -1
                        if {ra[2],ra[3]}&{r[2],r[3]}: return -1
                        c+=1
                    return c
                c1=cross(p1); c2=cross(p2)
                if c1!=1 or c2!=1: continue
                reach+=1
                # leftover tris of A not in a crossing must avoid B (checked in cross); now A's other half
                for a2 in PA[1-i]:
                    tb=tris(a2)
                    if tb&(t1|t2|ta) or KB in tb: continue
                    # open ends: triangle beyond each end unused
                    allt=t1|t2|ta|tb|{KA,KB}
                    ok=True
                    for p,K in ((p1,KB),(p2,KB),(a1,KA),(a2,KA)):
                        if beyond(p,K) in allt or not inH(beyond(p,K)): ok=False
                    if not ok: continue
                    hits+=1
                    for nm,p in (('Blo/Bhi',p1),('Blo/Bhi',p2),('A-crossing-half',a1),('A-free-half',a2)):
                        if p[2]: turnwhere[nm]+=1
                    break
    return hits,reach,turnwhere
if __name__=='__main__':
    L=int(sys.argv[1]) if len(sys.argv)>1 else 5
    for mode in ('strict','lesion'):
        h,r,tw=run(mode,L,6)
        print(mode,'L=',L,'realisations(B pair,A half,first A free half)=',h,' A-half crossing both B halves once=',r,' turns:',dict(tw))
