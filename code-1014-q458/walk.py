"""THE CENSUS IS A RANDOM WALK.

Claim:  L_{1-k}(z) = (z_1...z_k) * sum_{i=1}^k (z_i + 1/z_i).

Reason: in  sum over eps in {-1,0,1}^k with sum eps_i = 1-k, writing u=#{eps=+1},
v=#{eps=-1}, z=#{eps=0}: u-v = 1-k and u+v+z = k force z = 1-2u, hence u=0, z=1, v=k-1.
So EXACTLY ONE coordinate carries eps=0 (factor 1+z_i^2) and the rest carry eps=-1
(factor z_j).  Therefore

    L_{1-k}(z)^n = (z_1...z_k)^n * ( sum_i (z_i + 1/z_i) )^n

and  [x^{(1-k)(1^n)}] prod_i edot_{r_i}  =  W_n(r_1-n, ..., r_k-n),
the number of n-step walks on Z^k with steps +-e_i from 0 to (r-n).
"""
from itertools import product
from collections import defaultdict
from censusk import L_poly

def L_walkform(k):
    """(z_1..z_k) * sum_i (z_i+1/z_i), as dict of exponent tuples (may be negative)."""
    out = defaultdict(int)
    for i in range(k):
        for s in (+1,-1):
            key = tuple(1 + (s if j==i else 0) for j in range(k))
            out[key] += 1
    return dict(out)

def walk_counts(n, k, cap=None):
    """W_n(p) for all p, as dict; n-step walk on Z^k with steps +-e_i."""
    cur = {tuple([0]*k): 1}
    for _ in range(n):
        nxt = defaultdict(int)
        for p, c in cur.items():
            for i in range(k):
                for s in (+1,-1):
                    q = list(p); q[i] += s
                    nxt[tuple(q)] += c
        cur = dict(nxt)
    return cur

if __name__ == '__main__':
    print("=== L_{1-k} == (prod z)*sum(z_i+1/z_i) ? ===")
    for k in range(1,6):
        A = L_poly(k, 1-k)                       # built from the eps-sum definition
        B = L_walkform(k)                        # the closed form
        print(f"  k={k}: {A == B}   L = {sorted(A.items())[:4]}...")
    print("\n=== [z^r] L^n == W_n(r-n) ? ===")
    from censusk import poly_pow
    ok=fail=0
    for k in range(1,4):
        for n in range(1,7):
            Ln = poly_pow(L_poly(k,1-k), n, k, 2*n)
            W  = walk_counts(n,k)
            keys = set(Ln) | {tuple(p[i]+n for i in range(k)) for p in W}
            for r in keys:
                a = Ln.get(r,0); b = W.get(tuple(x-n for x in r),0)
                if a==b: ok+=1
                else: fail+=1; print("   MISMATCH",k,n,r,a,b)
    print(f"  ok={ok} fail={fail}")
