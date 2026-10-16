"""Q492: order of vanishing of Atilde_{M,I} on the long-root hyperplane
H_ab : w_a w_b = 1  (a<b),  as a function of the combinatorics of (M,I).

ord = multiplicity of the irreducible factor (w_a w_b - 1) in the numerator
      minus its multiplicity in the denominator.
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp, itertools
from collections import defaultdict
import q484_engine as E
q=E.q

def ord_H(f, w, a, b):
    """multiplicity of (w_a w_b - 1) in f (can be negative)."""
    f = sp.cancel(sp.together(f))
    if f == 0: return None
    num,den = sp.fraction(f)
    target = sp.expand(w[a-1]*w[b-1]-1)
    def mult(g):
        g = sp.expand(g)
        if g.is_number: return 0
        m=0
        for fac,k in sp.factor_list(g)[1]:
            fe=sp.expand(fac)
            if fe==target or fe==-target: m+=k
        return m
    return mult(num)-mult(den)

def main(N, X='C'):
    w=E.wvars(N); S=sorted(E.subsets(N), key=lambda M:(len(M),M))
    A={(M,I):E.Atilde(N,w,M,I,X) for M in S for I in S}
    print(f'=== N={N} type {X}: ord_{{w_a w_b=1}} Atilde_(M,I) ===')
    pairs=[(a,b) for a in range(1,N+1) for b in range(a+1,N+1)]
    recs=[]
    for (a,b) in pairs:
        for M in S:
            for I in S:
                if A[(M,I)]==0: continue
                o=ord_H(A[(M,I)],w,a,b)
                recs.append((a,b,M,I,o))
    # hypothesis: ord depends only on  |{a,b} cap M| and |{a,b} cap I|
    buckets=defaultdict(set)
    for a,b,M,I,o in recs:
        key=(len({a,b}&set(M)), len({a,b}&set(I)))
        buckets[key].add(o)
    print('\n  H1: ord determined by (|{a,b} n M|, |{a,b} n I|)?')
    ok=True
    for k in sorted(buckets):
        vals=sorted(buckets[k])
        print(f'    (|cap M|,|cap I|)={k}:  ord in {vals}')
        if len(vals)>1: ok=False
    print(f'  => {"DETERMINED" if ok else "NOT determined"} by H1')
    # refined: also record whether a,b in M as a *consecutive/order* datum
    print('\n  nonzero orders, listed:')
    for a,b,M,I,o in recs:
        if o!=0:
            print(f'    H_{a}{b}: M={M} I={I}  ord={o}')
    return recs

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 2)
