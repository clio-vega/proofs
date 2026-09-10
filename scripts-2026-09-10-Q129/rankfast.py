import sympy as sp
from ribbon import *
from repairB import compositions, g_comp
t=sp.Symbol('t')
def Zbeta(b):
    s=0;P=1
    for x in b: s+=x;P*=s
    return P
for n in (3,4):
    P=list(partitions(n)); C=list(compositions(n))
    A={}
    for b in C:
        for lam,c in g_comp(b,n).items(): A[(lam,b)]=sp.expand(c)
    rows=[];rhs=[]
    for i,lam in enumerate(P):
        for j,mu in enumerate(P):
            if j<i: continue
            rows.append([sp.expand(A.get((lam,b),0)*A.get((mu,b),0)) for b in C])
            rhs.append(1 if lam==mu else 0)
    M=sp.Matrix(rows); b=sp.Matrix(rhs)
    rM=M.rank(); rA=M.row_join(b).rank()
    M1=M.subs(t,-1); x1=sp.Matrix([sp.Rational(1,Zbeta(bb)) for bb in C])
    print(f"n={n} generic: rank(M)={rM} rank([M|b])={rA} -> {'INCONSISTENT' if rA>rM else 'consistent'}")
    print(f"      t=-1 CALIBRATION: known x=1/Z_beta satisfies system: {sp.simplify(M1*x1-b)==sp.zeros(len(rhs),1)}")
