"""Why do the SHORT roots (w_a = +-1) and the type-A roots (w_a = w_b) never
produce poles in c^nu, even though ord A_{nu,nu} > 0 there?

Measure ord A_{M,I} on those families and look for the rule.
Candidate (by the same admissible-K analysis):  ord_{w_a=1} A_{M,I} = [a in M][M\{a} not>= I]
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp
from collections import defaultdict
import q484_engine as E
q=E.q

def mult(g,t):
    g=sp.expand(g)
    if g.is_number or g==0: return 0
    m=0
    for fac,k in sp.factor_list(g)[1]:
        if sp.expand(fac-t)==0 or sp.expand(fac+t)==0: m+=k
    return m
def ordH(f,t):
    f=sp.cancel(sp.together(f))
    if f==0: return None
    n,d=sp.fraction(f); return mult(n,t)-mult(d,t)

for N in (2,3):
    w=E.wvars(N); S=sorted(E.subsets(N),key=lambda M:(len(M),M))
    A={(M,I):E.Atilde(N,w,M,I,'C') for M in S for I in S}
    print(f'=== N={N} ===')
    for sgn,lab in ((1,'w_a - 1'),(-1,'w_a + 1')):
        print(f'\n--- family {lab} ---  rule guess: [a in M][M\\{{a}} not>= I]')
        ok=True
        for a in range(1,N+1):
            T=sp.expand(w[a-1]-sgn)
            for M in S:
                for I in S:
                    if A[(M,I)]==0: continue
                    o=ordH(A[(M,I)],T)
                    Mm=tuple(x for x in M if x!=a)
                    pred = 1 if (a in M and not E.geq(Mm,I)) else 0
                    if o!=pred:
                        ok=False
                        print(f'   MISMATCH a={a} M={M} I={I}: measured {o}, guess {pred}')
        print(f'   guess {"CONFIRMED" if ok else "FAILED"}')
    print('\n--- family w_a - w_b (a<b) ---')
    for a in range(1,N+1):
        for b in range(a+1,N+1):
            T=sp.expand(w[a-1]-w[b-1])
            nz=[(M,I,ordH(A[(M,I)],T)) for M in S for I in S if A[(M,I)]!=0 and ordH(A[(M,I)],T)!=0]
            print(f'   H: w_{a}=w_{b}  nonzero ords: {[(M,I,o) for M,I,o in nz]}')
    print()
