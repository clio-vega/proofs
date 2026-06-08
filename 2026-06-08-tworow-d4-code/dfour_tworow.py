"""
dfour_tworow.py  --  shared machinery for the TWO-ROW d=4 order-law problem.

Background
----------
For the two-row shape lambda = (2m - b, b) the d=4 fibre value G_lambda(i) is purely
imaginary, and the order law for that family is equivalent to the non-vanishing of its
imaginary part I_b(m) for all integers m >= b.  Concretely

    G(u; m) = (1 - u)*(1 + s*u + u**2)**m,    s = 1 + i,
    I_b(m)  = Im( [u**b] G(u; m) )   in   ZZ[m].

I_b(m) has the forced rational roots m = 0, 1, ..., R with R = floor((b-1)/2):

    I_b(m) = (m)_{R+1} * Q_b(m),       (m)_{R+1} = m (m-1) ... (m-R),

deg_m Q_b = floor(b/2).  Writing q_b for the monic integer "irreducible part" of Q_b
(Q_b itself when 4 does not divide b; Q_b / (2m - (2b-1)) when 4 | b), the two-row d=4
law is equivalent to

        (diamond)   q_b has no rational root  (m >= b).

This module exposes:
    Ib(b)          -- I_b(m) as a sympy Poly in m over QQ (== ZZ scaled).
    Ib_direct(b,m0)-- I_b(m0) by an INDEPENDENT u-expansion (cross-check, no rl()).
    Qb(b)          -- Q_b(m) as a sympy Poly over QQ.
    qb_monic(b)    -- q_b, the monic primitive integer irreducible part, sympy Poly over ZZ.
    forced_factor(b)-- the product (m)_{R+1} as a Poly, for cross-checks.

Every builder is cross-checked against an independent evaluation before use.
"""

import sympy as sp
from sympy import symbols, expand, Poly, I, im, Rational

m = symbols('m', real=True)


# ----------------------------------------------------------------------------
#  I_b(m)  via the rl() expansion (the prior cycle's verified route)
# ----------------------------------------------------------------------------
def _rl(l):
    """Helper:  Im G = rl(b) - rl(b-1).  rl(l) = [u**l] of the m-th power piece."""
    s = 1 + I
    tot = 0
    for c in range(0, l // 2 + 1):
        k = l - 2 * c
        ff = 1
        for j in range(0, k + c):
            ff *= (m - j)
        tot += ff / (sp.factorial(k) * sp.factorial(c)) * s**k
    return expand(tot)


def Ib_expr(b):
    """I_b(m) as a sympy expression in m (integer-coefficient polynomial)."""
    return expand(im(_rl(b) - _rl(b - 1)))


def Ib(b):
    """I_b(m) as a Poly over QQ."""
    return Poly(Ib_expr(b), m, domain='QQ')


# ----------------------------------------------------------------------------
#  INDEPENDENT evaluation of I_b(m0):  expand (1-u)(1+s u+u^2)^{m0} directly.
#  Uses only integer m0 and a plain u-series -- shares no code with _rl().
# ----------------------------------------------------------------------------
def Ib_direct(b, m0):
    """Im [u**b] (1-u)(1+(1+i)u+u**2)**m0  for an explicit non-negative integer m0."""
    u = symbols('u')
    s = 1 + I
    series = expand((1 - u) * (1 + s * u + u**2)**m0)
    P = Poly(series, u)
    coeff = P.coeff_monomial(u**b) if b <= P.degree() else 0
    return int(im(expand(coeff)))


# ----------------------------------------------------------------------------
#  forced factor (m)_{R+1} and Q_b
# ----------------------------------------------------------------------------
def forced_factor(b):
    R = (b - 1) // 2
    P = Poly(1, m, domain='QQ')
    for r in range(0, R + 1):
        P = P * Poly(m - r, m, domain='QQ')
    return P


def Qb(b):
    """Q_b(m) = I_b(m) / (m)_{R+1}, as a Poly over QQ.  Exact division asserted."""
    P = Ib(b)
    R = (b - 1) // 2
    for r in range(0, R + 1):
        quo, rem = sp.div(P, Poly(m - r, m, domain='QQ'), m)
        assert rem.is_zero, f"forced root m={r} did not divide I_{b}"
        P = quo
    return P


def qb_monic(b):
    """Monic primitive integer irreducible part q_b of Q_b (Poly over ZZ)."""
    Q = Qb(b)
    facs = sp.factor_list(Q.as_expr())
    polys = [Poly(base, m) for base, _ in facs[1]]
    if not polys:
        return Poly(1, m, domain='ZZ')
    big = max(polys, key=lambda p: p.degree())
    return Poly(big.as_expr(), m, domain='ZZ')


# ----------------------------------------------------------------------------
#  self-test on import-as-main
# ----------------------------------------------------------------------------
def _selftest(bmax=16):
    print("cross-checking I_b via independent u-expansion ...")
    ok = True
    for b in range(1, bmax + 1):
        P = Ib(b)
        for m0 in (b, b + 1, b + 3, 2 * b + 1):
            lhs = int(P.eval(m0))
            rhs = Ib_direct(b, m0)
            if lhs != rhs:
                ok = False
                print(f"  MISMATCH b={b} m0={m0}: rl={lhs} direct={rhs}")
        # check forced factorisation reconstructs I_b
        recon = (forced_factor(b) * Qb(b)).as_expr()
        if expand(recon - Ib(b).as_expr()) != 0:
            ok = False
            print(f"  RECON FAIL b={b}")
    print("ALL CROSS-CHECKS PASSED" if ok else "FAILURES ABOVE")
    return ok


if __name__ == "__main__":
    _selftest()
