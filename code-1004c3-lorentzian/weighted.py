"""Why the dependency on CorollaryConvolution cannot be rerouted through the
doubly-corroborated labels (flow, CorollaryProduct, derivatives).

Those three operations act on NORMALISED coefficients.  Going through disjoint
variables and collapsing with `flow` yields prod_i N(Ptilde_i), whose raw
coefficient at gamma is the BINOMIALLY WEIGHTED count

    kw(gamma) = sum_{sum_i alpha^i = gamma}  gamma! / prod_i alpha^i!     (weights >= 1)

rather than the unweighted k(gamma) = #{...}.  Exhibit instances where the two
sequences differ, so that the two routes demonstrably prove different statements."""
import itertools
from fractions import Fraction
from lor import *

def Nof(c):
    import math
    return {a: Fraction(v, math.prod(math.factorial(t) for t in a)) for a,v in c.items()}
def pmulF(f,g):
    h={}
    for a,u in f.items():
        for b,v in g.items():
            k=tuple(x+y for x,y in zip(a,b)); h[k]=h.get(k,0)+u*v
    return h
import math
for (l,h) in [((0,0),(2,2)), ((1,0),(3,2)), ((0,0,0),(2,2,2))]:
    m=len(h); H=sum(h)
    f={(0,0,0):1}
    for i in range(m): f=poly_mult(f,Ptilde(l[i],h[i]))
    psi={(0,0,0):Fraction(1)}
    for i in range(m): psi=pmulF(psi,Nof(Ptilde(l[i],h[i])))
    # raw coefficients of the two routes, as sequences in a at fixed W-degree H-D
    for D in range(sum(l),H+1):
        k   = [f.get((a,D-a,H-D),0) for a in range(D+1)]
        kw  = [psi.get((a,D-a,H-D),Fraction(0))*math.factorial(a)*math.factorial(D-a)*math.factorial(H-D)
               for a in range(D+1)]
        kw  = [int(v) for v in kw]
        if k!=kw and any(k):
            print(f"  l={l} h={h} D={D}")
            print(f"    unweighted k (Theorem, route 2)       = {k}")
            print(f"    weighted  kw (flow+product route)     = {kw}")
            print(f"    both log-concave? k:{is_logconcave_pf2(k)}  kw:{is_logconcave_pf2(kw)}")
            break
print()
print("  Consistency check on the RECORDED form, settled by reasoning, not by re-reading:")
print("  if CorollaryConvolution were about f,g in DISJOINT variables then N(fg)=N(f)N(g)")
print("  identically, and the statement would be EXACTLY CorollaryProduct.  The record says")
print("  in terms that it is NOT an instance of CorollaryProduct.  So the recorded form is")
print("  internally coherent only under the SAME-variables reading -- which is the reading")
print("  Proposition stepB uses.  Verifying N(fg)=N(f)N(g) for disjoint variables:")
A=Ptilde(0,2)                      # in (X,E,W)
B=Ptilde(1,3)
def embed(c, lo, n):
    out={}
    for a,v in c.items():
        t=[0]*n
        for i,x in enumerate(a): t[lo+i]=x
        out[tuple(t)]=v
    return out
A6=embed(A,0,6); B6=embed(B,3,6)
lhs=Nof(poly_mult(A6,B6))
rhs=pmulF(Nof(A6),Nof(B6))
print("   N(fg) == N(f)N(g) on disjoint variables:", lhs==rhs, f"({len(lhs)} monomials)")
A3=Ptilde(0,2); B3=Ptilde(1,3)
lhs3=Nof(poly_mult(A3,B3)); rhs3=pmulF(Nof(A3),Nof(B3))
print("   N(fg) == N(f)N(g) on SHARED variables:  ", lhs3==rhs3, f"({len(lhs3)} monomials)")
