import sys, itertools, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B
t=sp.Symbol('t')
def sub(a,b):
    o={}
    for k in set(a)|set(b):
        x=sp.expand(a.get(k,0)-b.get(k,0))
        if x!=0: o[k]=x
    return o
def triple(lam,e,f,g):
    v={lam:sp.Integer(1)}; op=A.apply_op
    inner=lambda w: sub(op(op(w,g),f), op(op(w,f),g))
    return sub(op(inner(v),e), inner(op(v,e)))
print("claim: coefficient = -t^(e+g-5) (t^2-1)^2  on mu = unmaya(...),  e>=2,f>=3,g>=e+2, (e,g)!=(2,f+1)")
ok=bad=0
for e in range(2,5):
  for f in range(3,7):
    for g in range(e+2,9):
      if len({e,f,g})!=3: continue
      if e==2 and g==f+1: continue
      M=B.maya(()); b1,b2,b3=-e,-1,1-g
      if len({b1,b2,b3,b1+e,b2+f,b3+g})!=6: continue
      mu=B.unmaya(M-{b1,b2,b3}|{b1+e,b2+f,b3+g})
      pred=sp.expand(-t**(e+g-5)*(t**2-1)**2)
      act=triple((),e,f,g).get(mu,sp.Integer(0))
      good = sp.simplify(act-pred)==0
      ok+=good; bad+=not good
      if not good: print("  FAIL",(e,f,g),mu,"pred",pred,"act",act)
print(f"  {ok} confirmed, {bad} failed")
