"""THE LOAD-BEARING TEST.

Claim (twisted collapse on a long root).  For a<b with {a,b} subset M,
and H = H_ab = {w_a w_b = 1},

    Atilde_{M,I}  ==  rho * Atilde_{M\{a,b},I}   mod (w_a w_b - 1)

for EVERY I, with the UNIVERSAL multiplier

    rho = (1-q)^2 / [ (1 - q w_a^2)(1 - q w_a^{-2}) ] .

rho depends on neither M nor I nor K nor N.  Derived from Motegi's branching
formula by the observation that on w_b = w_a^{-1} every "one leg in {a,b}"
factor pairs with its own inverse and the ratio is 1.

Controls:
 (C1) PLANT A: a wrong multiplier rho*(1+q) must FAIL.
 (C2) PLANT B: the untwisted claim rho=1 must FAIL (that is thm:collision's
      shape, and the paper says long roots are NOT collision divisors).
 (C3) POSITIVE-DOMAIN check: report how many (M,I) pairs have BOTH sides
      nonzero -- a pair with both sides 0 cannot distinguish any multiplier.
"""
import sys; sys.path.insert(0,'.')
import sympy as sp, itertools, time
from collections import Counter
import q484_engine as E
q = E.q

def rho_universal(w,a,b):
    wa = w[a-1]
    return (1-q)**2 / ((1 - q*wa**2)*(1 - q/wa**2))

def divisible_by_H(expr, w, a, b):
    """is expr divisible by (w_a w_b - 1)?  i.e. does it vanish on H?"""
    expr = sp.cancel(sp.together(expr))
    if expr == 0: return True
    num,den = sp.fraction(expr)
    s = sp.Symbol('s')
    num_s = sp.expand(num.subs(w[b-1], s/w[a-1]))
    den_s = sp.expand(den.subs(w[b-1], s/w[a-1]))
    # den must not vanish at s=1 (else not a legitimate restriction)
    d1 = sp.simplify(den_s.subs(s,1))
    if d1 == 0: return None
    return sp.simplify(num_s.subs(s,1)) == 0

def run(N, X='C', mult=None, tag='CLAIM'):
    w = E.wvars(N); S = E.subsets(N)
    A = {(M,I): E.Atilde(N,w,M,I,X) for M in S for I in S}
    ok=0; bad=0; nondeg=0; undef=0; badlist=[]
    for a in range(1,N+1):
        for b in range(a+1,N+1):
            r0 = rho_universal(w,a,b) if mult is None else mult(w,a,b)
            for M in S:
                if a not in M or b not in M: continue
                Mp = tuple(x for x in M if x not in (a,b))
                for I in S:
                    lhs = A[(M,I)]; rhs = A[(Mp,I)]
                    d = divisible_by_H(lhs - r0*rhs, w, a, b)
                    # non-degenerate = at least one side nonzero on H
                    nz = (divisible_by_H(lhs,w,a,b) is False) or (divisible_by_H(rhs,w,a,b) is False)
                    if nz: nondeg += 1
                    if d is None: undef+=1
                    elif d: ok+=1
                    else:
                        bad+=1
                        if len(badlist)<6: badlist.append((a,b,M,Mp,I))
    print(f'  [{tag}] N={N} {X}: {ok} hold, {bad} FAIL, {undef} undefined'
          f'   (of which {nondeg} pairs are non-degenerate: some side != 0 on H)')
    for z in badlist: print('      fail at a,b,M,M\',I =',z)
    return ok,bad,nondeg

if __name__=='__main__':
    N = int(sys.argv[1]); X = sys.argv[2] if len(sys.argv)>2 else 'C'
    t0=time.time()
    run(N,X)
    run(N,X, mult=lambda w,a,b: rho_universal(w,a,b)*(1+q), tag='PLANT A (rho*(1+q))')
    run(N,X, mult=lambda w,a,b: sp.Integer(1),            tag='PLANT B (rho=1, untwisted)')
    print(f'  ({time.time()-t0:.0f}s)')
