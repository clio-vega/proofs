"""Q458 Task 1 -- the parameter census.

At k=1 (warn.py det size h=1), M=ell+2, N=2M:
    G = sum_{y in Z} t^{My^2-y} ( edot_{n-Ny} - edot_{n+2-Ny} )
    c_alpha(t) = [x^alpha] ( G * (x_1...x_n)^1 )
so for alpha=(1^n):
    c_{(1^n)}(t) = [x^0] G          (constant term, Laurent)

CLAIM (census of edot_r):
    [x^0] edot_r = binom(n, r/2) if r even and 0<=r<=2n, else 0.
because a monomial of e_r(x_1,1/x_1,...,x_n,1/x_n) is x^0 iff the chosen r-subset
is a union of |S| reciprocal PAIRS {x_i,1/x_i}, r=2|S|.

Hence with n=2m:
    c_{(1^n)}(t) = sum_{y in Z} t^{My^2-y} ( C(2m, m-My) - C(2m, m+1-My) ).
"""
import sys
from math import comb
import sympy as sp
t = sp.Symbol('t')

def closed_form(n, ell):
    """Predicted c_{(1^n)}(t) from the census."""
    M = ell + 2; N = 2*M
    if n % 2: return sp.Integer(0)
    m = n//2
    tot = sp.Integer(0)
    ymax = n + 4
    for y in range(-ymax, ymax+1):
        d = comb(n, m - M*y) if 0 <= m - M*y <= n else 0
        d -= comb(n, m+1 - M*y) if 0 <= m+1 - M*y <= n else 0
        if d: tot += d * t**(M*y*y - y)
    return sp.expand(tot)

def predicted_degree(n, ell):
    M = ell + 2; N = 2*M
    ys = [y for y in range(-(n+4), n+5) if -n <= N*y <= n+2]
    return max(M*y*y - y for y in ys), ys

if __name__ == '__main__':
    for ell in (0,1,2):
        for n in range(2, 13, 2):
            cf = closed_form(n, ell)
            d, ys = predicted_degree(n, ell)
            print(f"ell={ell} n={n} M={ell+2}  deg(closed form)={sp.Poly(cf,t).degree() if cf!=0 else None}"
                  f"  predicted_max={d}  admissible y={ys}")
