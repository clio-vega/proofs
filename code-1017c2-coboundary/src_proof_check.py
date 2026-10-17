"""Check the two intermediate ratios of the (SRC) proof separately -- a stronger test
than checking only that Atilde_s/Atilde_{s-1} = 1, because the two could both be wrong
and cancel.  Claim, type C (beta=1,gamma=0), eps^2=1, m in M, s=|M|, rest=M\{m}:

    F_s / F_{s-1} |_{w_m=eps}  =  prod_{k=1}^N (1-q eps w_k) * prod_{v in rest} (1-q eps/v)
    E_s / E_{s-1} |_{w_m=eps}  =  the SAME
  hence Atilde_{M,I}|_{w_m=eps} = Atilde_{M\{m},I}|_{w_m=eps}.
"""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1017c2-coboundary')
import sympy as sp, itertools
import q484_engine as E
q=E.q
def predicted(N,w,rest,eps):
    P=sp.Integer(1)
    for k in range(N): P*= (1-q*eps*w[k])
    for v in rest:     P*= (1-q*eps/v)
    return P
okF=okE=okA=0; badF=badE=badA=0; skip=0
for N in (2,3):
    w=E.wvars(N)
    for M in E.subsets(N):
        for m in M:
            Mm=tuple(j for j in M if j!=m)
            rest=[w[j-1] for j in Mm]
            for r in range(0,len(M)+1):
                for I in itertools.combinations(range(1,N+1),r):
                    for eps in (1,-1):
                        s=len(M)
                        if s-1 < r:   # F_{s-1} identically 0; ratio undefined, handled separately
                            skip+=1; continue
                        u1=[w[j-1] for j in M]; u0=[w[j-1] for j in Mm]
                        F1=E.F(N,s,r,u1,w,I,sp.Integer(1),sp.Integer(0)).subs(w[m-1],eps)
                        F0=E.F(N,s-1,r,u0,w,I,sp.Integer(1),sp.Integer(0)).subs(w[m-1],eps)
                        E1=E.EC(N,s,r,u1,w).subs(w[m-1],eps)
                        E0=E.EC(N,s-1,r,u0,w).subs(w[m-1],eps)
                        pred=predicted(N,w,rest,eps).subs(w[m-1],eps)
                        if sp.simplify(sp.cancel(E1/E0)-pred)==0: okE+=1
                        else: badE+=1; print('  E FAIL',N,M,m,I,eps)
                        if F0==0:
                            if sp.cancel(F1)==0: okF+=1
                            else: badF+=1; print('  F FAIL (F0=0,F1!=0)',N,M,m,I,eps)
                        elif sp.simplify(sp.cancel(F1/F0)-pred)==0: okF+=1
                        else: badF+=1; print('  F FAIL',N,M,m,I,eps)
                        A1=sp.cancel(E.Atilde(N,w,M,I,'C')).subs(w[m-1],eps)
                        A0=sp.cancel(E.Atilde(N,w,Mm,I,'C')).subs(w[m-1],eps)
                        if sp.simplify(sp.cancel(A1-A0))==0: okA+=1
                        else: badA+=1; print('  A FAIL',N,M,m,I,eps)
print(f'F-ratio  {okF} ok / {badF} bad')
print(f'E-ratio  {okE} ok / {badE} bad')
print(f'Atilde   {okA} ok / {badA} bad      (skipped s-1<r: {skip})')
# PLANT: the proof uses eps^2=1 in three places.  A non-root-of-unity eps must BREAK it.
print('--- plant: eps = 2 (not a root of unity) must break Atilde equality ---')
N=3; w=E.wvars(N); fired=0; tot=0
for M in [(1,2),(1,2,3),(2,3)]:
    for m in M:
        Mm=tuple(j for j in M if j!=m)
        for I in [(),(1,),(2,),(1,2)]:
            if len(Mm)<len(I): continue
            tot+=1
            A1=sp.cancel(E.Atilde(N,w,M,I,'C')).subs(w[m-1],2)
            A0=sp.cancel(E.Atilde(N,w,Mm,I,'C')).subs(w[m-1],2)
            if sp.simplify(sp.cancel(A1-A0))!=0: fired+=1
print(f'  plant fired {fired}/{tot} (must be >0; equality at eps=2 would mean the test is blind)')
