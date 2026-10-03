"""THEOREM M, tested in the ABSTRACT CLASS its hypotheses carve out -- not on slices.
(This is the 10-02 lesson: route (B) was 4,323,908/4,323,908 on slices and FALSE in its
own abstract class.)

Claim (M): for every box B = prod_i [P_i,Q_i] in Z^m, every sigma, every j,
   S = {y in B : sum_i y_i = sigma, sum_i (-y_i)_+ <= j}   is M-convex.
Also test the separable-convex generalisation with random convex phi_i, and the
NEGATIVE control: a NON-convex separable phi must break it."""
import mconv, itertools, random
from collections import Counter
random.seed(5)

def pts(B,sigma):
    return [y for y in itertools.product(*[range(P,Q+1) for (P,Q) in B]) if sum(y)==sigma]

def run(m,lo,hi,phi,label,maxbox=None):
    st=Counter(); ex=None
    rngs=[(P,Q) for P in range(lo,hi+1) for Q in range(P,hi+1)]
    boxes=list(itertools.product(rngs,repeat=m))
    if maxbox and len(boxes)>maxbox: boxes=random.sample(boxes,maxbox)
    for B in boxes:
        for sigma in range(sum(P for P,_ in B), sum(Q for _,Q in B)+1):
            Y=pts(B,sigma)
            if len(Y)<2: continue
            st['base']+=1
            if not mconv.exch_M(Y)[0]: st['CONTROL_BASE_BROKE']+=1
            vals=sorted({phi(y) for y in Y})
            for j in vals[:-1]:
                S=[y for y in Y if phi(y)<=j]
                if len(S)<2: continue
                st['tested']+=1; st['maxsize']=max(st['maxsize'],len(S))
                if mconv.exch_M(S)[0]: st['Mconvex']+=1
                else:
                    st['BROKE']+=1
                    if ex is None: ex=(m,B,sigma,j,sorted(S),mconv.exch_M(S)[1])
    print(f"  {label} m={m} box in [{lo},{hi}]: base sets {st['base']} "
          f"(base M-convex control broke {st['CONTROL_BASE_BROKE']}); sublevel tests {st['tested']}, "
          f"M-convex {st['Mconvex']}, BROKE {st['BROKE']}, max |S| {st['maxsize']}")
    if ex: print("     WITNESS:",ex)

kdef = lambda y: sum(max(-t,0) for t in y)
print("=== (M) with phi = k = sum_i (-y_i)_+  [the horizontal-strip defect] ===")
for m,lo,hi in [(2,-3,3),(3,-2,2),(3,-3,3),(4,-2,2)]:
    run(m,lo,hi,kdef,"k",maxbox=4000 if m>=4 else None)

print("=== (M) generalised: random SEPARABLE CONVEX phi ===")
for trial in range(4):
    slopes=[sorted(random.sample(range(-3,4),4)) for _ in range(3)]
    def phi(y,sl=slopes):
        tot=0
        for idx,t in enumerate(y):
            s=sl[idx%len(sl)]; v=0
            for q in range(min(t,0),0): v+= -s[0] if q< -1 else -s[1]
            for q in range(0,max(t,0)): v+= s[2] if q<1 else s[3]
            tot+=v
        return tot
    run(3,-2,2,phi,f"sepconvex#{trial}")

print("=== NEGATIVE control: NON-convex separable phi must BREAK (M) ===")
for trial in range(3):
    tab={t:random.randint(0,4) for t in range(-4,5)}
    nc=lambda y,tb=tab: sum(tb[t] for t in y)
    run(3,-2,2,nc,f"nonconvex#{trial}")
