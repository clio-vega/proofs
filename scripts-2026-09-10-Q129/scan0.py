import sympy as sp
from ribbon import *
t=sp.Symbol('t')
def preds(lam,N):
    out=[]
    for L in range(1,sum(lam)+1):
        for g,ht in remove_ribbons(lam,L,N+3): out.append((g,L,ht))
    return out

def system(n,mu):
    P=list(partitions(n)); pm={g:(L,ht) for g,L,ht in preds(mu,n)}
    gams=sorted(pm); rows=[];rhs=[]
    for lam in P:
        pl={g:(L,ht) for g,L,ht in preds(lam,n)}
        rows.append([t**pl[g][1] if (g in pl and pl[g][0]==pm[g][0]) else 0 for g in gams])
        rhs.append(1 if lam==mu else 0)
    return sp.Matrix(rows), sp.Matrix(rhs)

for n in (4,5,6,7):
    res={}
    for mu in partitions(n):
        M,b=system(n,mu)
        for val,nm in ((-1,'t=-1'),(0,'t=0'),(1,'t=1')):
            Ms,bs=M.subs(t,val),b.subs(t,val)
            ok = Ms.rank()==Ms.row_join(bs).rank()
            res.setdefault(nm,[]).append((mu,ok))
        ok=M.rank()==M.row_join(b).rank(); res.setdefault('generic',[]).append((mu,ok))
    print(f"n={n}:")
    for nm in ('t=-1','t=0','t=1','generic'):
        good=[m for m,o in res[nm] if o]; bad=[m for m,o in res[nm] if not o]
        print(f"   {nm:8s}: solvable for {len(good)}/{len(good)+len(bad)} mu" + (f"   FIRST FAILURE mu={bad[0]}" if bad else "   (all)"))
