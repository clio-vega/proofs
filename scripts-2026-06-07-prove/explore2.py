import sympy as sp
from sympy import symbols, factor, expand, Poly, I, im, Rational, ZZ, discriminant, primerange, gcd

m = symbols('m', real=True)

def ImP_coeff_b(b):
    s = 1 + I
    def rl(l):
        tot = 0
        for c in range(0, l//2 + 1):
            k = l - 2*c
            ff = 1
            for j in range(0, k + c):
                ff *= (m - j)
            tot += ff / (sp.factorial(k) * sp.factorial(c)) * s**k
        return sp.expand(tot)
    return expand(im(rl(b) - rl(b-1)))

def monic_part(b):
    """Return the monic integer polynomial q_b (the 'irreducible part': remove forced
       roots 0..floor((b-1)/2) and the half-integer linear factor when 4|b)."""
    Ib = ImP_coeff_b(b)
    P = Poly(Ib, m, domain='QQ')
    for r in range(0, (b-1)//2 + 1):
        P = sp.div(P, Poly(m - r, m, domain='QQ'), m)[0]
    # P is Q_b over QQ. Take its factor list, drop linear factors with non-integer root,
    # keep the monic integer poly of degree>=2 (there is exactly one nonlinear irreducible).
    facs = sp.factor_list(P.as_expr())
    pieces = []
    for base, mult in facs[1]:
        pb = Poly(base, m)
        pieces.append(pb)
    # The big factor is the one of max degree; make it monic-integer
    big = max(pieces, key=lambda p: p.degree())
    # clear to monic integer
    big = Poly(big.as_expr(), m, domain='QQ').monic()
    # convert to integer poly (should already be integer since it's a factor of integer-root poly)
    return Poly(big.as_expr(), m, domain='ZZ')

data = {}
for b in range(5, 31):
    q = monic_part(b)
    data[b] = q
    d = q.degree()
    c0 = q.all_coeffs()[-1]
    disc = sp.discriminant(q.as_expr(), m)
    # square-free part / is disc a perfect square?
    from sympy.ntheory.primetest import is_square as _issq
    is_perfect = _issq(abs(int(disc))) if disc != 0 else None
    print(f"b={b:2d} deg={d} const={c0} disc_is_square={is_perfect}")
