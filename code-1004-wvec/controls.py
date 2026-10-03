"""REAL controls for (H1''), with every prediction written down BEFORE the run.

(C0) remove one point from Y          -> must BREAK M-convexity  [tests the instrument]
(C1) SUPERlevel sets {k >= j}         -> predict BREAK (complement of a convex sublevel set)
(C2) phi(y)=sum_i |y_i - 1|  sublevel -> predict HOLD if the real theorem is
                                         'box cap hyperplane cap separable-convex sublevel'
(C3) phi(y)=#{i : y_i<0}     sublevel -> predict BREAK (phi not convex)
(C4) phi(y)=max_i |y_i|      sublevel -> convex but NOT separable; no prediction
Sizes reported: max and mean cardinality of the nontrivial sublevel sets tested."""
import gen, tent, mconv, random
from collections import Counter
def ys(mu,lam,n,m,b):
    t=sum(mu)+b; out=[]
    for nu in gen.box(mu,n,m):
        if sum(nu)!=t: continue
        lr=gen.LR(nu,lam,n,m)
        if any(R-L+1<=0 for (L,R) in lr): continue
        out.append(tuple(nu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)))
    return out
FUNS={'H1dd(k)': lambda y: sum(max(-t,0) for t in y),
      'C2(|y-1|)': lambda y: sum(abs(t-1) for t in y),
      'C3(#neg)' : lambda y: sum(1 for t in y if t<0),
      'C4(max)'  : lambda y: max(abs(t) for t in y)}
st=Counter(); sz=Counter(); ex={}
random.seed(11)
for m in (2,3,4,5):
  for n in range(m,m+6):
    for (mu,lam) in gen.pairs(n,m,9):
        d=sum(lam)-sum(mu)
        for b in range(0,d+1):
            Y=ys(mu,lam,n,m,b)
            if len(Y)<2: continue
            st['slices']+=1
            # C0: delete a point
            if len(Y)>=3:
                Z=list(Y); Z.pop(random.randrange(len(Z)))
                st['C0_broke' if not mconv.exch_M(Z)[0] else 'C0_still_Mconvex']+=1
            # C1: superlevel
            k=FUNS['H1dd(k)']; vals=sorted({k(y) for y in Y})
            for j in vals[1:]:
                S=[y for y in Y if k(y)>=j]
                if len(S)<2: continue
                st['C1_nontrivial']+=1
                st['C1_Mconvex' if mconv.exch_M(S)[0] else 'C1_broke']+=1
            for name,f in FUNS.items():
                vals=sorted({f(y) for y in Y})
                for j in vals[:-1]:
                    S=[y for y in Y if f(y)<=j]
                    if len(S)<2: continue
                    st[name+'_nontrivial']+=1
                    sz[name]=max(sz[name],len(S)); sz[name+'_tot']+=len(S)
                    if mconv.exch_M(S)[0]: st[name+'_Mconvex']+=1
                    else:
                        st[name+'_BROKE']+=1
                        ex.setdefault(name,(m,n,mu,lam,b,j,sorted(S),mconv.exch_M(S)[1]))
for k in sorted(st): print(f"  {k} = {st[k]}")
print("  sizes of nontrivial sublevel sets:",
      {n:(sz[n], round(sz[n+'_tot']/max(st[n+'_nontrivial'],1),2)) for n in FUNS}, "(max, mean)")
for k,v in sorted(ex.items()): print(" EX",k,v)
