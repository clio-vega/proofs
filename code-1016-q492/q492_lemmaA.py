"""LEMMA A test, all four root families.

Candidate hyperplanes are the non-monomial irreducible factors occurring in
ANY diagonal entry A_{I,I}.  For each, measure ord A_{M,I} for all (M,I) and
test the proposed rules.

Family (1) long root  H_ab : w_a w_b - 1   (a<b)   -- rule: [{a,b}<=M][M\{a,b} not>= I]
Families (2),(3),(4): report whether  ord A_{nu,nu} <= ord A_{nu,I}  for all I<=nu
(which would make the pole cancel automatically, with no condition on the jump).
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp
from collections import defaultdict
import q484_engine as E
q=E.q

def mult(g, target):
    g=sp.expand(g)
    if g.is_number or g==0: return 0
    m=0
    for fac,k in sp.factor_list(g)[1]:
        fe=sp.expand(fac)
        if sp.expand(fe-target)==0 or sp.expand(fe+target)==0: m+=k
    return m

def ordH(f,target):
    f=sp.cancel(sp.together(f))
    if f==0: return None
    n,d=sp.fraction(f)
    return mult(n,target)-mult(d,target)

def run(N):
    w=E.wvars(N); S=sorted(E.subsets(N),key=lambda M:(len(M),M))
    A={(M,I):E.Atilde(N,w,M,I,'C') for M in S for I in S}
    # collect candidate hyperplanes from all diagonal entries
    cand={}
    for I in S:
        f=sp.cancel(sp.together(A[(I,I)]))
        for part in sp.fraction(f):
            for fac,k in sp.factor_list(sp.expand(part))[1]:
                if len(sp.Poly(fac,*w).monoms())>1:
                    cand[sp.srepr(sp.expand(fac))]=sp.expand(fac)
    print(f'=== N={N}: {len(cand)} candidate hyperplanes from diagonal entries ===')
    def classify(t):
        s=sp.sstr(sp.factor(t))
        return s
    # long-root family targets
    longroot={}
    for a in range(1,N+1):
        for b in range(a+1,N+1):
            longroot[(a,b)]=sp.expand(w[a-1]*w[b-1]-1)
    print('\n--- LEMMA A on the long-root family  H_ab: w_a w_b = 1 ---')
    ok=True; checked=0
    for (a,b),T in longroot.items():
        for M in S:
            for I in S:
                if A[(M,I)]==0: continue
                o=ordH(A[(M,I)],T); checked+=1
                Mm=tuple(x for x in M if x not in (a,b))
                pred = 1 if (a in M and b in M and not E.geq(Mm,I)) else 0
                if o!=pred:
                    ok=False
                    print(f'   MISMATCH H_{a}{b} M={M} I={I}: measured {o}, rule {pred}')
    print(f'   rule [{{a,b}}<=M][M\\{{a,b}} not>= I] : {"CONFIRMED" if ok else "FAILED"}  ({checked} pairs (M,I) tested)')
    # other families
    print('\n--- other candidate hyperplanes: is  ord A_{nu,nu} <= ord A_{nu,I}  for all I with nu>=I ? ---')
    for key,T in sorted(cand.items(), key=lambda kv: sp.sstr(kv[1])):
        if any(sp.expand(T-t)==0 for t in longroot.values()): continue
        bad=[]
        for nu in S:
            dnn=ordH(A[(nu,nu)],T)
            if dnn is None or dnn<=0: continue
            for I in S:
                if A[(nu,I)]==0: continue
                d=ordH(A[(nu,I)],T)
                if d<dnn: bad.append((nu,I,dnn,d))
        tag = 'OK (auto-cancel)' if not bad else f'FAILS ({len(bad)} cases, e.g. {bad[0]})'
        print(f'   {sp.sstr(sp.factor(T)):26s}  {tag}')
    return

for N in (2,3):
    run(N); print()
