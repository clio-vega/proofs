"""STEP 0.  Is the BASIS itself regular on the A-type divisors {w_i = w_j}?

Atilde_{M,I} is built from F (a sum over S_n) which contains factors
(q u_j - u_k)/(u_j - u_k).  With u = w_M, if i,j in M then some u_j - u_k IS
w_i - w_j, so EVERY SINGLE SUMMAND has an A-type pole.  The question is whether
the symmetrization kills it.  Measure, do not assume.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1018-atype')
import sympy as sp, q484_engine as E, divisors as D
q = E.q

def main(N, X='C'):
    w = E.wvars(N); S = sorted(E.subsets(N), key=lambda M:(len(M),M))
    t0=time.time()
    tally = {}
    polar = []
    for M in S:
        for I in S:
            v = E.Atilde(N, w, M, I, X)
            if v == 0: continue
            kinds = D.den_classify(v, w, q)
            for k, lst in kinds.items():
                tally[k] = tally.get(k,0)+1
                if k == 'A':
                    polar.append((M,I,lst))
    print(f'=== N={N} type {X}: basis entries Atilde_(M,I), {time.time()-t0:.0f}s ===')
    print(f'    nonzero entries: {sum(1 for M in S for I in S if E.Atilde(N,w,M,I,X)!=0)}')
    print(f'    denominator buckets over nonzero entries: {tally if tally else "ALL LAURENT (no non-monomial denominators at all)"}')
    print(f'    A-TYPE POLES: {len(polar)}')
    for M,I,lst in polar: print(f'      M={M} I={I}: {lst}')

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 3, sys.argv[2] if len(sys.argv)>2 else 'C')
