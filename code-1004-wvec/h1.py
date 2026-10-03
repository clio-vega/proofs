"""(H1): is W = {w(nu) : nu in effective slice} M-convex?
Structural obstruction known in advance: every M-convex set lies in ONE hyperplane
sum_i w_i = const (exchange w-e_i+e_j preserves the coordinate sum), while
A-general-m-tent gives sum_i w_i(nu) = 2*Lambda(nu)+m.  So (H1) fails on every slice
carrying two distinct half-widths.  Measure how often that is, and check the rest."""
import gen, tent, mconv, sys
from collections import Counter

def Wset(mu,lam,n,m,b):
    target=sum(mu)+b; out=[]
    for nu in gen.box(mu,n,m):
        if sum(nu)!=target: continue
        lr=gen.LR(nu,lam,n,m); w=tuple(R-L+1 for (L,R) in lr)
        if any(x<=0 for x in w): continue
        out.append(w)
    return out

PLAN=[(2,range(2,11),10),(3,range(3,10),9),(4,range(4,9),8),(5,range(5,9),8)]
st=Counter(); ex={}
for m,nrange,dmax in PLAN:
  for n in nrange:
    G=n-m
    for (mu,lam) in gen.pairs(n,m,dmax):
        d=sum(lam)-sum(mu)
        for b in range(0,d+1):
            ws=Wset(mu,lam,n,m,b)
            if not ws: continue
            S=sorted(set(ws)); st['slices']+=1
            hw=sorted({(sum(w)-m)//2 for w in S}); nsum=len(hw)
            st['Mprof_len_%d'%min(nsum,6)]+=1; st['G=%d'%G]+=1; st['m=%d'%m]+=1
            st['nW_%d'%min(len(S),8)]+=1
            if nsum>1:
                st['H1_FALSE_by_sum']+=1
                if 'sum' not in ex: ex['sum']=(m,n,mu,lam,b,S,mconv.sums(S))
            else:
                ok,wit=mconv.exch_M(S)
                st['oneplane_Mconvex' if ok else 'oneplane_NOT_Mconvex']+=1
                if not ok and 'oneplane_bad' not in ex: ex['oneplane_bad']=(m,n,mu,lam,b,S,wit)
print("slices:",st['slices'])
print(" (H1) FALSE by the coordinate-sum obstruction alone:",st['H1_FALSE_by_sum'],
      "= %.2f%%"%(100*st['H1_FALSE_by_sum']/st['slices']))
print(" single-hyperplane slices: M-convex",st['oneplane_Mconvex']," NOT",st['oneplane_NOT_Mconvex'])
print(" M-profile length:",{k:v for k,v in sorted(st.items()) if k.startswith('Mprof_len')})
print(" |W| (set size):",{k:v for k,v in sorted(st.items()) if k.startswith('nW_')})
print(" G:",{k:v for k,v in sorted(st.items()) if k.startswith('G=')})
print(" m:",{k:v for k,v in sorted(st.items()) if k.startswith('m=')})
for k,v in ex.items(): print(" EX",k,v)
