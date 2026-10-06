"""Corollary 4 of math/0202090, tested exhaustively, with the non-vacuity gate.

    I_alpha(u,w)  =  sum_v  c^w_{u,v} * I_alpha(w_0 v, w_0)

Conventions fixed and VALIDATED (see validate_thm2): labeled Bruhat order
u --(k,b)--> w on covers with u^{-1}w = (i,j), i <= k < j, b = u(i) = w(j);
chains strictly increasing in lex on (k,b); type = exponent vector of first
coordinates; w_0 v = w_0 . v (left multiplication), the Poincare dual
convention c^{w_0}_{u, w_0 u} = 1.
"""
import sys, collections, itertools
from ls import *
from schub import schubert, structure_constant, structure_constants_linalg, pmul, length as _l

def cor4_data(n):
    W0=w0(n)
    # I-tables for the right-hand side pieces, once
    Itab_rhs={}
    for v in perms(n):
        Itab_rhs[v]=I_table(mul(W0,v), W0, n)
    # structure constants
    C=collections.defaultdict(dict)   # C[(u,w)][v] = c^w_{u,v}
    for u in perms(n):
        for v in perms(n):
            m=length(u)+length(v)
            if m>length(W0): continue
            for w in perms(n):
                if length(w)!=m: continue
                c=structure_constant(u,v,w,n)
                if c: C[(u,w)][v]=c
    return W0, Itab_rhs, C

def test_cor4(n, verbose=True):
    W0,Itab_rhs,C = cor4_data(n)
    checks=0; fails=[]
    falsifiable=0; triples=0
    inv_fail=[]
    per_triple=[]
    for u in perms(n):
        for w in perms(n):
            m=length(w)-length(u)
            if m<0: continue
            lhs=I_table(u,w,n)
            cs=C.get((u,w),{})
            # ---- invariant check: total cardinality, per alpha and overall
            tot_l=sum(len(x) for x in lhs.values())
            tot_r=sum(c*sum(len(x) for x in Itab_rhs[v].values()) for v,c in cs.items())
            if tot_l!=tot_r: inv_fail.append((u,w,tot_l,tot_r))
            alphas=set(lhs)|{a for v in cs for a in Itab_rhs[v]}
            triples+=len(alphas)
            for al in alphas:
                L=len(lhs.get(al,[]))
                R=sum(c*len(Itab_rhs[v].get(al,[])) for v,c in cs.items())
                checks+=1
                if L!=R: fails.append((u,w,al,L,R))
                # NON-VACUITY: |Gamma_alpha(u,w)|>1 AND >=2 distinct v with c>0
                nv = [v for v,c in cs.items() if c>0 and len(Itab_rhs[v].get(al,[]))>0]
                if L>1 and len(nv)>=2:
                    falsifiable+=1
                    per_triple.append((u,w,al,L,{v:cs[v] for v in nv}))
    if verbose:
        print('n=%d'%n)
        print('  (u,w,alpha) instances of Corollary 4 checked : %d' % checks)
        print('  disagreements                                : %d' % len(fails))
        for f in fails[:5]: print('     FAIL u=%s w=%s alpha=%s lhs=%d rhs=%d'%f)
        print('  total-cardinality invariant failures         : %d' % len(inv_fail))
        for f in inv_fail[:5]: print('     INV u=%s w=%s |Gamma|=%d rhs=%d'%f)
        print('  *** FALSIFIABLE instances (|Gamma_alpha(u,w)|>1 AND >=2 v) : %d of %d'
              % (falsifiable, checks))
    return dict(checks=checks,fails=fails,falsifiable=falsifiable,
                inv_fail=inv_fail,per_triple=per_triple,C=C,Itab_rhs=Itab_rhs)

if __name__=='__main__':
    for n in (2,3,4):
        test_cor4(n); print()
