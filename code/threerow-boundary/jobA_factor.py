"""
jobA_factor.py -- Job A step 3: hunt for the structural 2^k factor of master N_i.

Factor the EVEN- and ODD-slice polynomials N_i(P,B) over Q (these fit reliably by
tensor interpolation, unlike the parity-constrained (a,b) fit).  On the slice,
   a = 2P + b + c,  a-b+1 = 2P+c+1,  a-c+2 = 2P+2B+2 = 2(P+B+1),  b = 2B (or 2B+1).
For each irreducible factor report its 2-adic slice content (min v2 over a grid).
We want to see WHERE the fixed-divisor content comes from: an even-spaced run of
linear factors, or an explicit 2^k*(odd) block -- the sibling of Theorem B's
(a-b+1) factor.
"""
import sympy as sp
from math import factorial
from fast_alt import Malt

P, B = sp.symbols('P B')


def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def Ni_val(c, i, Pv, Bv, even):
    k = c - i
    bv = 2 * Bv if even else 2 * Bv + 1
    av = 2 * Pv + bv + c
    m = Pv + bv + c
    j = bv + i
    if bv < c + 1 or j > m:
        return None
    M = Malt(av, bv, c, j, m)
    den = (bv - c + 1)
    for t in range(2, i + 1):
        den *= (bv + t)
    num = M * factorial(c) * factorial(k)
    if den == 0 or num % den != 0:
        return None
    return num // den


def fit_slice_poly(c, i, even, pdeg, bdeg):
    monos = [(p, q) for p in range(pdeg + 1) for q in range(bdeg + 1)]
    B0 = ((c + 1) // 2) + 2
    pts, rhs = [], []
    for Pv in range(0, pdeg + 1):
        for Bv in range(B0, B0 + bdeg + 1):
            v = Ni_val(c, i, Pv, Bv, even)
            if v is None:
                return None
            pts.append((Pv, Bv))
            rhs.append(sp.Integer(v))
    A = sp.Matrix([[sp.Integer(Pv)**p * sp.Integer(Bv)**q for (p, q) in monos]
                   for (Pv, Bv) in pts])
    sol = A.LUsolve(sp.Matrix(rhs))
    return sp.expand(sum(sol[t] * P**p * B**q for t, (p, q) in enumerate(monos)))


def factor_content(fac, c, even, grid=96):
    g = 10**9
    B0 = ((c + 1) // 2) + 1
    fe = sp.lambdify((P, B), sp.expand(fac), modules='sympy')
    for Pv in range(grid):
        for Bv in range(B0, grid):
            val = int(fe(Pv, Bv))
            if val == 0:
                continue
            vv = v2int(val)
            if vv < g:
                g = vv
    return g


def degrees_for(k):
    return (k // 2 + 3, k + 4)


CASES = [(6, 2, 4), (6, 3, 3), (5, 2, 3), (7, 2, 5), (8, 2, 6),
         (7, 4, 3), (9, 4, 5)]

if __name__ == "__main__":
    for (c, i, k) in CASES:
        print("=" * 84)
        print("c=%d i=%d k=%d   master N_i(P,B);  a=2P+b+c, a-b+1=2P+c+1, "
              "a-c+2=2(P+B+1)" % (c, i, k))
        pdeg, bdeg = degrees_for(k)
        for even in (True, False):
            expr = fit_slice_poly(c, i, even, pdeg, bdeg)
            slc = "EVEN (b=2B)" if even else "ODD (b=2B+1)"
            cst, facs = sp.factor_list(expr)
            num = sp.Rational(cst)
            cst_v2 = v2int(num.p) - v2int(num.q)
            print("  [%s]  leading const = %s  (v2=%d)" % (slc, cst, cst_v2))
            for fac, mult in facs:
                fc = factor_content(fac, c, even)
                print("     x%d  %-44s  content=%d  (total %d)"
                      % (mult, sp.factor(fac), fc, fc * mult))
