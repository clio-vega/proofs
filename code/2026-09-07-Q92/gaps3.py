import sys, itertools, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B, engineC as C
t=sp.Symbol('t')
def assignments(M,Mp,e,f,g):
    R=sorted(M-Mp); Ad=sorted(Mp-M); out=[]
    for pr in itertools.permutations(R):
        b1,b2,b3=pr
        if sorted([b1+e,b2+f,b3+g])==Ad and len({b1,b2,b3,b1+e,b2+f,b3+g})==6:
            out.append(pr)
    return out
rows=[]
for e in range(2,6):
  for f in range(3,7):
    for g in range(e+2,9):
      if len({e,f,g})!=3: continue
      M=B.maya(()); b1,b2,b3=-e,-1,1-g
      if len({b1,b2,b3,b1+e,b2+f,b3+g})!=6: continue
      Mp=M-{b1,b2,b3}|{b1+e,b2+f,b3+g}
      n=len(assignments(M,Mp,e,f,g))
      rows.append((e,f,g,n))
uniq=[r for r in rows if r[3]==1]; multi=[r for r in rows if r[3]!=1]
print(f"unique assignment: {len(uniq)}/{len(rows)}")
print("NON-unique cases:", [(e,f,g,n) for e,f,g,n in multi])
