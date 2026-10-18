"""STEP 6.  Test Theorem A's CONCLUSION directly, not just its two inputs.

Theorem A says  L_ij = G_ij.  I have checked both colengths (V4, orbit count), but
that is the PROOF's arithmetic.  Here I test the equality itself, both inclusions,
as a two-arm test:

  ARM G->L (the nontrivial inclusion).  Take f in G_ij and expand it in the phi
  basis: f = sum_I a_I phi_I.  Theorem A says every a_I lies in O_ij, i.e.
  ord_(w_i-w_j) a_I >= 0.  Generators of G_ij used: the s_ij-symmetric indicators
  1_{M} + 1_{sM}, the uniformiser times a single indicator t.1_{M}, and 1_{M} for
  s_ij-FIXED M.  These generate G_ij over O_ij.
  CONTROL (must FIRE): 1_{M} for a MOVED M is NOT in G_ij (it is 1 at M, 0 at sM,
  and 1 - 0 = 1 is not in tO).  Some a_I must then have ord = -1.  If this does not
  fire, the test cannot see the difference between L and G and the green arm is
  decoration.

  ARM L->G is Lemmas A1+A2 and is already checked in V2/V3.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1018-atype')
import sympy as sp
from itertools import combinations
import q484_engine as E
q = E.q

def ordA(f, w, i, j):
    f = sp.cancel(sp.together(f))
    if f == 0: return sp.oo
    num, den = sp.fraction(f)
    def mult(p):
        k=0; p=sp.expand(p)
        while True:
            qq, rr = sp.div(sp.Poly(p, w[j-1]), sp.Poly(w[j-1]-w[i-1], w[j-1]))
            if sp.expand(rr.as_expr())!=0: return k
            p=qq.as_expr(); k+=1
    return mult(num)-mult(den)

def swap(M,i,j):
    s=set(M)
    if (i in s)!=(j in s): s.symmetric_difference_update({i,j})
    return tuple(sorted(s))

def expand(vec, S, A, w):
    """solve sum_I a_I Atilde_{M,I} = vec(M) by triangularity (same solve as the
    structure constants, so a bug there would show up here too)."""
    a={}
    for nu in S:
        acc=sp.Integer(0)
        for kap in S:
            if kap==nu: break
            if a.get(kap,0)!=0: acc+=a[kap]*A[(nu,kap)]
        a[nu]=sp.cancel(sp.together((vec[nu]-acc)/A[(nu,nu)]))
    return a

def main(N=3, X='C'):
    w=E.wvars(N); S=sorted(E.subsets(N), key=lambda M:(len(M),M))
    A={(M,I):sp.cancel(E.Atilde(N,w,M,I,X)) for M in S for I in S}
    print(f'=== N={N} type {X}: Theorem A conclusion, L_ij =?= G_ij ===', flush=True)
    for i,j in combinations(range(1,N+1),2):
        t = w[i-1]-w[j-1]
        moved=[M for M in S if len({i,j}&set(M))==1]
        fixed=[M for M in S if len({i,j}&set(M))!=1]
        gens=[]
        # orbit sums (in G)
        seen=set()
        for M in moved:
            Ms=swap(M,i,j)
            if M in seen: continue
            seen|={M,Ms}
            gens.append(('orbit-sum %s+%s'%(M,Ms), {K:(sp.Integer(1) if K in (M,Ms) else sp.Integer(0)) for K in S}, True))
            gens.append(('t.1_%s'%(M,),          {K:(t if K==M else sp.Integer(0)) for K in S}, True))
            gens.append(('CONTROL 1_%s (moved)'%(M,), {K:(sp.Integer(1) if K==M else sp.Integer(0)) for K in S}, False))
        for M in fixed[:2]:
            gens.append(('1_%s (s-fixed)'%(M,), {K:(sp.Integer(1) if K==M else sp.Integer(0)) for K in S}, True))
        okc=badc=0; ctl_fire=ctl_tot=0
        for name, vec, in_G in gens:
            a=expand(vec,S,A,w)
            ords=[ordA(a[I],w,i,j) for I in S]
            neg=[o for o in ords if o!=sp.oo and o<0]
            if in_G:
                if neg: badc+=1; print(f'   *** ({i},{j}) {name}: in G but coefficient ord {neg} < 0 -- THEOREM A FALSE', flush=True)
                else: okc+=1
            else:
                ctl_tot+=1
                if neg: ctl_fire+=1
                else: print(f'   ({i},{j}) {name}: NOT in G yet all coefficients regular -- control did not fire', flush=True)
        print(f'  (i,j)=({i},{j}): G-generators with all coefficients in O: {okc}/{okc+badc}'
              f'{"  *** %d VIOLATIONS ***"%badc if badc else ""}'
              f'   |  CONTROL (1_M for moved M, NOT in G): fires {ctl_fire}/{ctl_tot}'
              f' -> {"FIRES" if ctl_fire==ctl_tot and ctl_tot else "BLIND"}', flush=True)

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 3, sys.argv[2] if len(sys.argv)>2 else 'C')
