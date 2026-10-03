"""The repair the diagnosis points to.

Diagnosis: the half-width defect k(y) = sum_i (-y_i)_+ is a CONVEX function of y
(A-half-width-is-l1-distance).  The strata the brief's (H1)/(H1') use are its LEVEL
sets -- shells of a convex function -- which is why they have holes, and the shells sit
at even distance in sum_i w_i, which is the parity obstruction.
Its SUBLEVEL sets are the convex objects.  And they live in y-space, where sum_i y_i = sigma
is CONSTANT, so M-convexity is well posed.

(H1''): for every j >= 0, Y_{<=j} = {y in Y : k(y) <= j} is M-convex.
Controls: Y_{<=infty} = Y = box cap hyperplane must come back M-convex (positive control);
a perturbed stratification (k' = sum_i (y_i)_+, and k''= sum_i |y_i|) is the refusal test."""
import gen, tent, mconv
from collections import Counter

def ys(mu,lam,n,m,b):
    target=sum(mu)+b; out=[]
    for nu in gen.box(mu,n,m):
        if sum(nu)!=target: continue
        lr=gen.LR(nu,lam,n,m)
        if any(R-L+1<=0 for (L,R) in lr): continue
        out.append(tuple(nu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)))
    return out

def kdef(y):  return sum(max(-t,0) for t in y)
def kalt(y):  return sum(max(t,0) for t in y)
def kabs(y):  return sum(abs(t) for t in y)

PLAN=[(2,range(2,12),12),(3,range(3,11),10),(4,range(4,10),9),(5,range(5,9),8),(6,range(6,9),7)]
st=Counter(); ex={}
for m,nrange,dmax in PLAN:
  for n in nrange:
    for (mu,lam) in gen.pairs(n,m,dmax):
        d=sum(lam)-sum(mu)
        for b in range(0,d+1):
            Y=ys(mu,lam,n,m,b)
            if not Y: continue
            st['slices']+=1
            if len({sum(y) for y in Y})!=1: st['SIGMA_VARIES']+=1
            # positive control: the whole Y (box cap hyperplane) must be M-convex
            if len(Y)>=2:
                st['ctrl_nontrivial']+=1
                st['ctrl_Y_Mconvex' if mconv.exch_M(Y)[0] else 'ctrl_Y_NOT']+=1
                if not mconv.exch_M(Y)[0]: ex.setdefault('ctrl',(m,n,mu,lam,b,sorted(Y)))
            for name,f in [('H1dd',kdef),('refuse_plus',kalt),('refuse_abs',kabs)]:
                K=sorted({f(y) for y in Y})
                for j in K[:-1]:          # j = max is the whole Y, already the control
                    S=sorted([y for y in Y if f(y)<=j])
                    if len(S)<2: st[name+'_trivial']+=1; continue
                    st[name+'_nontrivial']+=1
                    if mconv.exch_M(S)[0]: st[name+'_Mconvex']+=1
                    else:
                        st[name+'_NOT']+=1
                        ex.setdefault(name,(m,n,mu,lam,b,j,S,mconv.exch_M(S)[1]))
for k in sorted(st): print(f"  {k} = {st[k]}")
print()
for k,v in sorted(ex.items()): print(" EX",k,v)
