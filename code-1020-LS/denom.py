r"""V5b: the load-bearing half of Lemma 2.1(2):
     den( Atilde^X_{A,B} ) |_{q=0}  is nonzero and involves NO w_c with c not in A.
(The lemma does not claim it is a monomial -- only that it is q-free and
A-supported, which is exactly what Lemma 2.1(3),(4) consume.)
DOMAIN: N = 1..3 exhaustive, both types, every (A,B) with A >= B.
"""
import sympy as sp, sys
sys.path.insert(0,'.'); import q484_engine as E
q = E.q
bad=[]; ok=0; mono=0
for N in range(1,4):
    w = E.wvars(N)
    for X in ('C','B'):
        for A in E.subsets(N):
            for B in E.subsets(N):
                if not E.geq(A,B): continue
                v = sp.cancel(sp.together(E.Atilde(N,w,A,B,X=X)))
                num,den = sp.fraction(v)
                d0 = sp.expand(den.subs(q,0))
                if d0 == 0:
                    bad.append(('den|q=0 VANISHES',N,X,A,B)); continue
                p = sp.Poly(d0, *w)
                outside = [j+1 for j in range(N)
                           if (j+1) not in A and any(mm[j]!=0 for mm in p.monoms())]
                if outside:
                    bad.append((N,X,A,B,sp.factor(d0),'uses',outside))
                else:
                    ok += 1
                    if len(p.monoms())==1: mono += 1
print(f"den|_{{q=0}} nonzero and supported on (w_a)_{{a in A}}: {ok} ok, {len(bad)} bad")
print(f"   (of the {ok}, {mono} are in fact single monomials -- stronger than the lemma claims)")
for b in bad[:6]: print("  ",b)
