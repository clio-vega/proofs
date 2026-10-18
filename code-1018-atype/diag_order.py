"""STEP 2.  The colength arithmetic.

CLAIM (derived, then tested).  Fix i<j.  By (M4) the diagonal Atilde_{I,I} is the
root product diag_rootproduct, whose factors are
   (1 - 1/(w_a w_b)), a<=b in I            pair / boundary divisors
   (1 - w_a/w_{i_l}),  a in I^c, i_l in I, i_l > a      <-- the only A-TYPE factors
   (1 - q w_a w_b), (1 - q w_a/w_{i_l})     q-deformed
A factor (1 - w_a/w_{i_l}) vanishes on {w_i = w_j} iff {a, i_l} = {i,j}; with i<j and
the constraint i_l > a this forces a = i, i_l = j.  Hence
   ord_{(w_i-w_j)} Atilde_{I,I}  =  1 if (i not in I and j in I), else 0,
and summing over the 2^N labels
   sum_I ord = #{I : i notin I, j in I} = 2^{N-2}.
That is EXACTLY the colength of the A-type GKM lattice in O^P (one condition per
2-element s_ij-orbit, and there are 2^{N-2} of them).  Tested here directly against
the engine, not against the root product -- an INDEPENDENT mechanism (S_n sum vs
product over positive roots), which is the control that caught the E^B error.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1018-atype')
import sympy as sp
from itertools import combinations
import q484_engine as E
q = E.q

def ord_at(f, w, i, j):
    """ord_{(w_i - w_j)}(f) for f a rational function, via repeated substitution."""
    f = sp.cancel(sp.together(f))
    if f == 0: return sp.oo
    num, den = sp.fraction(f)
    def mult(p):
        k = 0
        p = sp.expand(p)
        while True:
            qq, rr = sp.div(sp.Poly(p, w[j-1]), sp.Poly(w[j-1]-w[i-1], w[j-1]))
            if sp.expand(rr.as_expr()) != 0: return k
            p = qq.as_expr(); k += 1
    return mult(num) - mult(den)

def main(X='C', Ns=(2,3,4)):
    for N in Ns:
        w = E.wvars(N); S = sorted(E.subsets(N), key=lambda M:(len(M),M))
        t0=time.time()
        print(f'=== N={N} type {X} ===', flush=True)
        # cross-check (M4): engine diagonal == root product
        agree = all(sp.cancel(E.Atilde(N,w,I,I,X) - E.diag_rootproduct(N,w,I))==0 for I in S) if X=='C' else 'N/A (no type-B root product in the engine; ord measured directly below)'
        print(f'  (M4) engine diagonal == product over positive roots, all {len(S)} labels: {agree}', flush=True)
        for i,j in combinations(range(1,N+1),2):
            tot = 0; rows=[]
            for I in S:
                o = ord_at(E.Atilde(N,w,I,I,X), w, i, j)
                pred = 1 if (i not in I and j in I) else 0
                tot += o
                if o != pred: rows.append((I,o,pred))
            print(f'  (i,j)=({i},{j}): sum_I ord = {tot}, predicted 2^(N-2) = {2**(N-2)}'
                  f'  {"OK" if tot==2**(N-2) else "*** MISMATCH ***"}'
                  f'   per-label formula ord=[i notin I, j in I]: {"all agree" if not rows else "*** %s ***"%rows}', flush=True)
        print(f'  ({time.time()-t0:.0f}s)', flush=True)

if __name__=='__main__':
    main(sys.argv[1] if len(sys.argv)>1 else 'C',
         tuple(int(x) for x in sys.argv[2:]) or (2,3,4))
