"""General-k census -- Task 2 / Q457.

For alpha=(1^n), c_alpha = [x^{(1^n)}](G * prodx^k) = [x^{b.(1,..,1)}] G with b = 1-k.

GENERALISED CENSUS.  For a product of k edots and a CONSTANT target exponent vector
b.(1,...,1), the n variables decouple:

    [x^{b(1..1)}] prod_{i=1}^k edot_{r_i}  =  [z_1^{r_1}...z_k^{r_k}]  L_b(z_1..z_k)^n
    L_b(z) = sum_{eps in {-1,0,1}^k, sum eps_i = b}  prod_i ( z_i if eps_i!=0 else 1+z_i^2 )

because for each variable x_j and each factor i the local choice is a subset of
{x_j, 1/x_j}: empty (exp 0, size 0), {x_j} (+1,1), {1/x_j} (-1,1), both (0,2).
k=1,b=0 gives L=1+z^2 -> binom(n,r/2), recovering the k=1 census.
"""
from itertools import product, permutations
from collections import defaultdict

def L_poly(k, b):
    """L_b as dict {(d_1..d_k): coeff}."""
    out = defaultdict(int)
    for eps in product((-1,0,1), repeat=k):
        if sum(eps) != b: continue
        # expand prod_i (z_i if eps_i else 1+z_i^2)
        terms = [[(1,1)] if e else [(0,1),(2,1)] for e in eps]
        for pick in product(*terms):
            key = tuple(p[0] for p in pick)
            out[key] += 1
    return dict(out)

def poly_pow(P, n, kvars, cap):
    """P^n as dict, truncating any exponent > cap (cap = 2n, max edot index)."""
    R = {tuple([0]*kvars): 1}
    base, e = dict(P), n
    while e:
        if e & 1: R = poly_mul(R, base, cap)
        e >>= 1
        if e: base = poly_mul(base, base, cap)
    return R

def poly_mul(A, B, cap):
    out = defaultdict(int)
    for a, ca in A.items():
        for b_, cb in B.items():
            key = tuple(x+y for x, y in zip(a, b_))
            if any(v > cap for v in key): continue
            out[key] += ca*cb
    return {k2: v for k2, v in out.items() if v}

def c_1n(n, k, ell, verbose=False):
    """c_{(1^n)}(t) as dict {t-exponent: int}, general k.  No sympy, no determinant
    expansion in x -- the x-side is done entirely by the census."""
    M = k + ell + 1; N = 2*M; b = 1 - k
    Ln = poly_pow(L_poly(k, b), n, k, 2*n)
    def CT(rs):
        if any(r < 0 or r > 2*n for r in rs): return 0
        return Ln.get(tuple(rs), 0)
    ymax = (2*n + 2*k + 4)//N + 1
    out = defaultdict(int)
    for ys in product(range(-ymax, ymax+1), repeat=k):
        for sigma in permutations(range(1, k+1)):
            sgn = perm_sign(sigma)
            texp = sum(M*ys[i]*ys[i] - sigma[i]*ys[i] for i in range(k))
            # prod_i ( edot_{n-(i+1)+sigma_i-N y_i} - edot_{n+(i+1)+sigma_i-N y_i} )
            tot = 0
            for S in product((0,1), repeat=k):
                rs = [n + (1 if S[i] else -1)*(i+1) + sigma[i] - N*ys[i] for i in range(k)]
                tot += (-1)**sum(S) * CT(rs)
            if tot: out[texp] += sgn*tot
    return {e: v for e, v in out.items() if v}

def perm_sign(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i+1:
            j = p.index(i+1); p[i], p[j] = p[j], p[i]; s = -s
    return s
