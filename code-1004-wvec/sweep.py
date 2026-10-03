"""THE SWEEP.  Four well-posed questions about the width-vector set W of a real slice.

Q0 (Lemma P): is sum_i w_i = G+m-sum|y_i|, and is its parity constant = (G+m+sigma) mod 2?
Q1 (H1, the brief): is W M-convex?            -- predicted FALSE whenever >=2 half-widths
Q1nat: is W M^natural-convex?                 -- predicted FALSE for the same reason
Q2 (H1', the level repair): is each level W_k = {w : k(nu)=k} M-convex?  (constant sum there)
Q3 (the graded lift): is hat W = {(k(nu), w(nu))} M^nat-convex?  (unit-step sums, parity gone)

Object sizes reported throughout (hazard 5): |W|, number of levels, max level size,
diameter of W, G = n-m, m.  Every count conditioned on the test being NON-TRIVIAL."""
import gen, tent, mconv, sys
from collections import Counter, defaultdict

def rows(mu,lam,n,m,b):
    target=sum(mu)+b; out=[]
    for nu in gen.box(mu,n,m):
        if sum(nu)!=target: continue
        lr=gen.LR(nu,lam,n,m); w=tuple(R-L+1 for (L,R) in lr)
        if any(x<=0 for x in w): continue
        y=tuple(nu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m))
        out.append((y,w))
    return out

PLAN=[(2,range(2,12),12),(3,range(3,11),10),(4,range(4,10),9),(5,range(5,9),8),(6,range(6,9),7)]
st=Counter(); ex={}
for m,nrange,dmax in PLAN:
  for n in nrange:
    G=n-m
    for (mu,lam) in gen.pairs(n,m,dmax):
        d=sum(lam)-sum(mu)
        for b in range(0,d+1):
            R=rows(mu,lam,n,m,b)
            if not R: continue
            st['slices']+=1; st['m=%d'%m]+=1; st['G=%d'%G]+=1
            sig=sum(R[0][0])
            # ---- Q0
            pars=set()
            for y,w in R:
                if sum(w)==G+m-sum(abs(t) for t in y): st['Q0_id_ok']+=1
                else:
                    st['Q0_id_BAD']+=1; ex.setdefault('Q0id',(m,n,mu,lam,b,y,w))
                pars.add(sum(w)%2)
            if pars=={(G+m+sig)%2}: st['Q0_parity_ok']+=1
            else: st['Q0_parity_BAD']+=1; ex.setdefault('Q0par',(m,n,mu,lam,b,pars,(G+m+sig)%2))
            W=sorted({w for y,w in R})
            c=max((G-sig),0)/2
            lev=defaultdict(set)
            for y,w in R: lev[sum(max(-t,0) for t in y)].add(w)
            st['|W|=%s'%min(len(W),9)]+=1
            st['nlev=%s'%min(len(lev),6)]+=1
            st['maxlev=%s'%min(max(len(v) for v in lev.values()),6)]+=1
            # ---- Q1 / Q1nat : only meaningful when |W|>=2
            if len(W)>=2:
                st['Q1_nontrivial']+=1
                if len(mconv.sums(W))>1: st['Q1_multisum']+=1
                if mconv.exch_M(W)[0]: st['Q1_Mconvex']+=1; ex.setdefault('Q1yes',(m,n,mu,lam,b,W))
                else: st['Q1_not']+=1
                if mconv.exch_Mnat(W)[0]: st['Q1nat_yes']+=1; ex.setdefault('Q1natyes',(m,n,mu,lam,b,W))
                else: st['Q1nat_not']+=1
            # ---- Q2 : levels; only meaningful for levels of size >=2
            for k,S in lev.items():
                S=sorted(S)
                if len(S)<2: st['Q2_level_trivial']+=1; continue
                st['Q2_level_nontrivial']+=1
                if len(mconv.sums(S))!=1:
                    st['Q2_level_SUMVARIES']+=1; ex.setdefault('Q2sum',(m,n,mu,lam,b,k,S))
                if mconv.exch_M(S)[0]: st['Q2_level_Mconvex']+=1
                else:
                    st['Q2_level_NOT']+=1
                    ex.setdefault('Q2no',(m,n,mu,lam,b,k,S,mconv.exch_M(S)[1]))
            # ---- Q3 : graded lift
            hatW=sorted({(k,)+w for k,S in lev.items() for w in S})
            if len(hatW)>=2:
                st['Q3_nontrivial']+=1
                s=mconv.sums(hatW)
                if s==list(range(s[0],s[-1]+1)): st['Q3_sums_interval']+=1
                else: st['Q3_sums_GAPPED']+=1; ex.setdefault('Q3gap',(m,n,mu,lam,b,hatW,s))
                if mconv.exch_Mnat(hatW)[0]: st['Q3_Mnat_yes']+=1
                else:
                    st['Q3_Mnat_not']+=1
                    ex.setdefault('Q3no',(m,n,mu,lam,b,hatW,mconv.exch_Mnat(hatW)[1]))
for k in sorted(st): print(f"  {k} = {st[k]}")
print()
for k,v in sorted(ex.items()): print(" EX",k,v)
