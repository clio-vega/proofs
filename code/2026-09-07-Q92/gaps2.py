import sys, itertools, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B, engineC as C
t=sp.Symbol('t')
def nm(lam,mu): return len(B.maya(lam)^B.maya(mu))//2
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

print("=== uniform bead-3 family: lambda=empty, b1=-e (moves e), b2=-1 (f), b3=1-g (g)")
print("    valid for e>=2, f>=3, g>=e+2 ; predicted coefficient  -(t - 1/t)^2 * t^P")
ok=bad=0
for e in range(2,5):
  for f in range(3,6):
    for g in range(e+2,8):
      if len({e,f,g})!=3: continue
      M=B.maya(())
      b1,b2,b3=-e,-1,1-g
      sites={b1,b2,b3,b1+e,b2+f,b3+g}
      if len(sites)!=6: print("  site collision",e,f,g); bad+=1; continue
      if not all(x in M for x in (b1,b2,b3)) or any(x in M for x in (b1+e,b2+f,b3+g)):
          print("  illegal",e,f,g); bad+=1; continue
      Mp=M-{b1,b2,b3}|{b1+e,b2+f,b3+g}
      mu=B.unmaya(Mp)
      P=C.cnt(M,b1,b1+e)+C.cnt(M,b2,b2+f)+C.cnt(M,b3,b3+g)
      def kap(bi,ei,bj,ej): return (1 if bj<bi+ei<bj+ej else 0)-(1 if bj<bi<bj+ej else 0)
      k32=kap(b3,g,b2,f); k21=kap(b2,f,b1,e); k31=kap(b3,g,b1,e); K=k21+k31
      pred=sp.expand(t**P*(t**k32-t**(-k32))*(t**K-t**(-K)))
      act=triple((),e,f,g).get(mu,sp.Integer(0))
      good = (k32!=0 and K!=0 and nm((),mu)==3 and sp.simplify(act-pred)==0 and sp.expand(act)!=0)
      ok+=good; bad+=(not good)
      if not good: print("  FAIL",(e,f,g),mu,"k32",k32,"K",K,"pred",pred,"act",act)
print(f"    {ok} confirmed nonzero bead-3 elements, {bad} failed")
