"""THEOREM B (Q457).  Completing the square in the t-exponent.

With q_i := N y_i and p_i := sigma(i) + eps_i i - N y_i  (eps_i = +1 iff i in S),

    E(y,sigma) = M sum y_i^2 - sum sigma(i) y_i
               = (1/4M) [ sum (q_i - sigma(i))^2 - sum sigma(i)^2 ]
               = (1/4M) [ sum_i (p_i - eps_i i)^2 - sum_{s=1}^k s^2 ].

sigma CANCELS (it survives only in the constant sum s^2).  So maximising the
t-exponent over all contributing configurations is maximising
    sum_i (|p_i| + i)^2     subject to    ||p||_1 <= n
(the walk-support constraint), a CONVEX max over an l1 ball: unique vertex
p = +- n e_k, value (n+k)^2 + sum_{s<k} s^2, giving

    E_max = n(n+2k)/(4M),   realisable iff N | n  (eps = -1^k, y = -(n/N) e_k, coeff +1)
                                       or N | n+2k (eps_k=+1, y = ((n+2k)/N) e_k, coeff -1).
"""
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1014-q458')
from fractions import Fraction
from censusk import c_1n

def predict(n,k,ell):
    M=k+ell+1; N=2*M
    if n % 2: return None
    if n % N == 0:      return n*(n+2*k)//(4*M), +1, 'B'
    if (n+2*k) % N == 0: return n*(n+2*k)//(4*M), -1, 'A'
    return None

if __name__=='__main__':
    ok=fail=skip=0
    rows=[]
    for k in range(1,6):
        for ell in range(0,4):
            M=k+ell+1; N=2*M
            nmax = {1:36,2:30,3:24,4:18,5:14}[k]
            for n in range(2,nmax+1,2):
                pr = predict(n,k,ell)
                if pr is None: skip+=1; continue
                D,sign,br = pr
                d = c_1n(n,k,ell)
                if not d: print(f"  !! k={k} ell={ell} n={n}: c==0"); fail+=1; continue
                gotD=max(d); gotc=d[gotD]
                good = (gotD==D and gotc==sign)
                ok+= good; fail+= (not good)
                if not good: print(f"  MISMATCH k={k} ell={ell} n={n} branch={br}: pred deg={D} coef={sign}, got deg={gotD} coef={gotc}")
                rows.append((k,ell,n,br,D,sign,gotD,gotc,good))
    print(f"THEOREM B: ok={ok} fail={fail}  (skipped {skip} n outside the two residue classes)")
    # non-vacuity / satisfiability control: print a sample of the numbers themselves
    print("\nsample rows (k,ell,n,branch,predicted deg,predicted coef,measured deg,measured coef):")
    for r in rows[::max(1,len(rows)//25)]: print("  ",r)
