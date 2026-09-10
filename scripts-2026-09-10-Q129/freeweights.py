import sympy as sp
from ribbon import *
from repairB import compositions, g_comp
t = sp.Symbol('t')
def Zbeta(b):
    s=0;P=1
    for x in b: s+=x;P*=s
    return P

print("=== T4a: the omega-cocycle reciprocity for the composition matrix ===")
print("    claim  A(t)_{lam',beta} == t^{n-l(beta)} A(1/t)_{lam,beta}")
for n in (2,3,4,5):
    P=list(partitions(n)); C=list(compositions(n))
    A={}
    for b in C:
        for lam,c in g_comp(b,n).items(): A[(lam,b)]=c
    ok=True
    for b in C:
        for lam in P:
            l=sp.expand(A.get((conj(lam),b),0))
            r=sp.expand(t**(n-len(b))*sp.together(A.get((lam,b),0)).subs(t,1/t))
            if sp.simplify(l-r)!=0: ok=False; print("   FAIL",n,lam,b,l,r)
    print(f"  n={n}: {ok}")

print()
print("=== T4b: FREE the weights.  Does ANY diagonal x_beta(t) give A(t) diag(x) A(t)^T = I ? ===")
for n in (3,4,5):
    P=list(partitions(n)); C=list(compositions(n))
    A={}
    for b in C:
        for lam,c in g_comp(b,n).items(): A[(lam,b)]=sp.expand(c)
    xs=sp.symbols(f'x0:{len(C)}')
    for target,tname in ((None,'I'),('P','omega')):
        eqs=[]
        for i,lam in enumerate(P):
            for j,mu in enumerate(P):
                if j<i: continue
                lhs=sum(xs[k]*A.get((lam,b),0)*A.get((mu,b),0) for k,b in enumerate(C))
                tgt = (1 if lam==mu else 0) if target is None else (1 if conj(lam)==mu else 0)
                eqs.append(sp.expand(lhs-tgt))
        sol=sp.solve(eqs, xs, dict=True)
        print(f"  n={n} #unknowns={len(C)} #equations={len(eqs)}  target={tname}: solution over Q(t) -> {sol if sol else 'NONE'}")
