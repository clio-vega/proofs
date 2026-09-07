import sys, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B
t=sp.Symbol('t')
def nm(lam,mu): return len(B.maya(lam)^B.maya(mu))//2

print("=== TEST 1: is C_{e,f}=[R_e,R_f]/(1+t) a MULTIPLICATION operator? ===")
for (e,f) in [(1,2),(1,3),(2,3)]:
    c0={k:sp.cancel(v/(1+t)) for k,v in A.commutator((),e,f).items()}
    c1={k:sp.cancel(v/(1+t)) for k,v in A.commutator((1,),e,f).items()}
    # if C = M_g then g = C(s_empty); predict C(s_1) = g * s_1
    g=c0                      # g as a Schur expansion
    pred={}
    for lam,co in g.items():
        for mu,w in A.R_e(lam,1).items():   # Pieri: g*s_1 = sum over adding a box
            pred[mu]=sp.expand(pred.get(mu,0)+co*w)
    pred={k:v for k,v in pred.items() if sp.expand(v)!=0}
    same=set(pred)==set(c1) and all(sp.simplify(pred[k]-c1[k])==0 for k in c1)
    print(f"  (e,f)=({e},{f})  g=C(s_0)={g}")
    print(f"      C(s_1) actual = {c1}")
    print(f"      g*s_1        = {pred}")
    print(f"      MULTIPLICATION OPERATOR? {same}")

print("\n=== TEST 2: is [R_e,R_f] in span{R_g}? ===")
for (e,f) in [(1,2),(1,3),(2,3)]:
    com=A.commutator((),e,f); Rg={g:A.R_e((),g) for g in range(1,e+f+1)}
    print(f"  ({e},{f}): [R_e,R_f]s_0 supported on {sorted(com)}")
    print(f"           R_{e+f} s_0 supported on {sorted(Rg[e+f])}")

print("\n=== TEST 3: does the body-number stay bounded?  [R_2,[R_3,R_4]] on s_0 ===")
def applyop(vec,e): return A.apply_op(vec,e)
v={():sp.Integer(1)}
inner={}   # [R_3,R_4] s_0
a=applyop(applyop(v,4),3); b=applyop(applyop(v,3),4)
for k in set(a)|set(b):
    x=sp.expand(a.get(k,0)-b.get(k,0))
    if x!=0: inner[k]=x
outer={}
a2={};b2={}
for lam,co in inner.items():
    for mu,w in A.R_e(lam,2).items(): a2[mu]=sp.expand(a2.get(mu,0)+co*w)
c2=applyop(inner,2)
d2={}
tmp=applyop(v,2)   # R_2 s_0
# [R_2, X] s_0 = R_2 (X s_0) - X (R_2 s_0)
xr={}
aa=applyop(applyop(tmp,4),3); bb=applyop(applyop(tmp,3),4)
for k in set(aa)|set(bb):
    x=sp.expand(aa.get(k,0)-bb.get(k,0))
    if x!=0: xr[k]=x
for k in set(c2)|set(xr):
    x=sp.expand(c2.get(k,0)-xr.get(k,0))
    if x!=0: outer[k]=x
bodies={}
for mu,co in outer.items(): bodies.setdefault(nm((),mu),[]).append((mu,co))
for r in sorted(bodies): print(f"   {r}-bead sector: {len(bodies[r])} terms, e.g. {bodies[r][0]}")
