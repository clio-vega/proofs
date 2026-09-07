import sys, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B, engineC as C
ok=bad=0; fails=[]
NMAX=int(sys.argv[1]); EMAX=int(sys.argv[2])
for e in range(1,EMAX+1):
  for f in range(1,EMAX+1):
    for n in range(0,NMAX+1):
      for lam in A.partitions(n):
        a=A.commutator(lam,e,f); c=C.routes(lam,e,f)
        keys=set(a)|set(c)
        if all(sp.simplify(sp.expand(a.get(k,0)-c.get(k,0)))==0 for k in keys): ok+=1
        else: bad+=1; fails.append((e,f,lam,a,c))
print(f"A vs C (route count):  AGREE {ok}/{ok+bad}   |lam|<={NMAX}, 1<=e,f<={EMAX}  [includes e=f]")
for x in fails[:5]: print("FAIL",x)
