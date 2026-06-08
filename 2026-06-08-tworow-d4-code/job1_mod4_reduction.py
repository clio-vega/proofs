"""
job1_mod4_reduction.py  --  JOB 1.1 of the d=4 two-row CODE plan.

For b = 0 (mod 4), verify EXACTLY (as a congruence of integer-valued polynomials
in m, i.e. all Mahler / forward-difference coefficients agree mod 4) that

    I_b(m) == -m*C(m-1,R) + 2*C(m,2)*C(m-2,R) - 2*C(m,3)*C(m-3,R-1)   (mod 4),

with R = floor((b-1)/2) = b/2 - 1.

WHY this is the right test.  Two integer-valued polynomials f, g satisfy
f(m) == g(m) (mod 4) for every integer m  IFF  their Mahler expansions
f = sum c_k C(m,k), g = sum d_k C(m,k) have c_k == d_k (mod 4) for all k.
The Mahler coefficients are the forward differences c_k = (Delta^k f)(0), a FINITE
integer set (k = 0..deg).  So checking them mod 4 is a rigorous, finite proof of the
polynomial congruence -- not a sample of values.

We also expose the MECHANISM: I_b is a signed sum over j of C(m;j,c)*Im((1+i)^j)
(c determined by j+2c in {b, b-1}); we confirm that every j >= 5 (odd) and j >= 6
(even) carries Im((1+i)^j) == 0 (mod 4), so only j in {1,2,3} survive -- which is
exactly the displayed three-term reduction.
"""

import sympy as sp
from sympy import symbols, binomial, expand, Integer
from dfour_tworow import Ib, m

def mahler_coeffs(expr):
    """Forward differences c_k = (Delta^k f)(0) for k = 0..deg(expr) in m.
    expr must be an integer-valued polynomial in m; returns list of python ints."""
    P = sp.Poly(expr, m)
    d = P.degree()
    vals = [Integer(P.eval(i)) for i in range(d + 1)]
    coeffs = []
    for k in range(d + 1):
        # k-th forward difference at 0
        ck = sum((-1)**(k - i) * sp.binomial(k, i) * vals[i] for i in range(k + 1))
        coeffs.append(int(ck))
    return coeffs

def rhs_reduction(b):
    """The displayed three-term mod-4 reduction, as a sympy expression."""
    R = (b - 1) // 2
    e = (-m * binomial(m - 1, R)
         + 2 * binomial(m, 2) * binomial(m - 2, R)
         - 2 * binomial(m, 3) * binomial(m - 3, R - 1))
    return sp.expand(sp.expand_func(e))  # turn binomials into polynomials in m

def im_s_pow(j):
    """Im((1+i)^j) as an integer."""
    val = (1 + sp.I)**j
    return int(sp.im(sp.expand(val)))

def main():
    print("=" * 78)
    print("JOB 1.1  --  exactness of the mod-4 reduction for b = 0 (mod 4)")
    print("=" * 78)
    all_ok = True
    for b in range(8, 49, 4):
        R = (b - 1) // 2
        D = expand(Ib(b).as_expr() - rhs_reduction(b))
        cks = mahler_coeffs(D)
        bad = [(k, c) for k, c in enumerate(cks) if c % 4 != 0]
        status = "EXACT (mod 4)" if not bad else f"FAILS at Mahler coeffs {bad}"
        if bad:
            all_ok = False
        print(f"  b={b:2d}  R={R:2d}  deg(diff)={len(cks)-1:2d}  -> {status}")
    print()
    print("-" * 78)
    print("MECHANISM:  Im((1+i)^j) mod 4 for the surviving/​killed j")
    print("-" * 78)
    print("  j :  Im((1+i)^j)   (mod 4)   note")
    for j in range(0, 13):
        v = im_s_pow(j)
        note = ""
        if j in (1, 2, 3):
            note = "<- SURVIVES (in reduction)"
        elif j == 0 or j % 4 == 0:
            note = "Im = 0 (j = 0 mod 4)"
        elif v % 4 == 0:
            note = "killed mod 4"
        print(f"  {j:2d} : {v:6d}        {v % 4:2d}      {note}")
    print()
    print("Surviving j in {1,2,3} give exactly:")
    print("   j=1 (from -u part): -C(m;1,R)        = -m*C(m-1,R)")
    print("   j=2 (from  1 part): +2*C(m;2,R)      = +2*C(m,2)*C(m-2,R)")
    print("   j=3 (from -u part): -2*C(m;3,R-1)    = -2*C(m,3)*C(m-3,R-1)")
    print()
    print("RESULT:", "ALL b in {8,...,48} EXACT mod 4." if all_ok else "SOME FAILED -- see above.")
    return all_ok

if __name__ == "__main__":
    main()
