"""
Computational fingerprinting of I_b(m) and its deflated part Q_b(m).

I_b(m) = Im( [u^b]( (1-u)*(1 + s*u + u^2)^m ) ),  s = 1+i.
       = C_b(m) - C_{b-1}(m),  C_l(m) = Im([u^l](1+su+u^2)^m).

C_l(m) = sum_{c>=0, a=l-2c>=0} [ falling(m, a+c)/(a! c!) ] * Im((1+i)^a)

I_b(m) = (prod_{r=0}^{floor((b-1)/2)} (m-r)) * Q_b(m).
"""
import sympy as sp
from sympy import symbols, I, factorial, factorint, expand, Rational, Poly, lcm, gcd, prod
from functools import reduce

m = symbols('m', real=True)


def falling(x, k):
    """falling factorial x*(x-1)*...*(x-k+1), k factors."""
    r = sp.Integer(1)
    for j in range(k):
        r *= (x - j)
    return r


def C_l(l):
    """Im([u^l](1+s u+u^2)^m) as exact poly in m, s=1+i."""
    total = sp.Integer(0)
    c = 0
    while l - 2 * c >= 0:
        a = l - 2 * c
        imval = sp.im((1 + I) ** a)
        if imval != 0:
            term = falling(m, a + c) / (factorial(a) * factorial(c)) * imval
            total += term
        c += 1
    return sp.expand(total)


def I_b(b):
    return sp.expand(C_l(b) - C_l(b - 1))


def primitive_int_poly(poly_expr, var=m):
    """Return (content_rational, primitive_poly_expr) where primitive has integer
    coprime coeffs and positive leading coeff."""
    p = Poly(sp.expand(poly_expr), var)
    coeffs = p.all_coeffs()
    if all(c == 0 for c in coeffs):
        return sp.Integer(0), sp.Integer(0)
    # clear denominators
    denoms = [sp.nsimplify(c).as_numer_denom()[1] for c in coeffs]
    L = reduce(lcm, denoms, sp.Integer(1))
    int_coeffs = [sp.Integer(c * L) for c in coeffs]
    g = reduce(gcd, [abs(c) for c in int_coeffs if c != 0], sp.Integer(0))
    int_coeffs = [c / g for c in int_coeffs]
    # positive leading
    if int_coeffs[0] < 0:
        int_coeffs = [-c for c in int_coeffs]
    prim = sum(int_coeffs[k] * var ** (len(int_coeffs) - 1 - k)
               for k in range(len(int_coeffs)))
    return int_coeffs, sp.expand(prim)


def deflate(b):
    """Divide out forced roots {0..floor((b-1)/2)}, verify exact, return Q_b expr."""
    Ib = I_b(b)
    forced = list(range(0, (b - 1) // 2 + 1))
    Q = Ib
    for r in forced:
        # verify root
        val = Ib.subs(m, r)
        q, rem = sp.div(sp.Poly(Q, m), sp.Poly(m - r, m))
        Q = q.as_expr()
    # verify Q * prod = Ib
    recon = Q * prod([(m - r) for r in forced])
    assert sp.expand(recon - Ib) == 0, f"deflation mismatch b={b}"
    return Ib, Q, forced


if __name__ == "__main__":
    pass
