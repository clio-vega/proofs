import sympy as sp
from ribbon import *
from repairB import compositions, g_comp
t = sp.Symbol('t')

def Zbeta(b):
    s=0; P=1
    for x in b: s+=x; P*=s
    return P

for n in (2,3,4,5):
    P=list(partitions(n)); idx={p:i for i,p in enumerate(P)}
    C=list(compositions(n))
    A={}
    for b in C:
        v=g_comp(b,n)
        for lam,c in v.items(): A[(lam,b)]=c
    M=sp.zeros(len(P),len(P))
    for lam in P:
        for mu in P:
            s=0
            for b in C:
                a1=A.get((lam,b),0); a2=A.get((mu,b),0)
                if a1!=0 and a2!=0: s+= sp.Rational(1,Zbeta(b))*a1*a2
            M[idx[lam],idx[mu]]=sp.expand(s)
    Mm1=M.subs(t,-1)
    print(f"n={n}: A(t)Z^-1A(t)^T == I  at t=-1 ? {Mm1==sp.eye(len(P))}")
    gen = sp.simplify(M-sp.eye(len(P)))==sp.zeros(len(P),len(P))
    print(f"      ... for general t ? {gen}")
    if not gen and n<=4:
        for lam in P:
            for mu in P:
                v=sp.factor(M[idx[lam],idx[mu]]-(1 if lam==mu else 0))
                if v!=0: print(f"        ({lam},{mu}): defect = {v}")
