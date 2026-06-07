import sympy as sp
from sympy import symbols, factor, expand, Poly, I, im, Rational, ZZ, QQ, discriminant, primerange

m = symbols('m', real=True)

def ImP_coeff_b(b):
    """ I_b(m) = [u^b]((1-u) Im P), Im P = sum_{k odd}(-1)^((k-1)/2) C(m,k) u^k W^{m-k}.
        Build as polynomial in m using the rl_sym style (Im G = rl(b)-rl(b-1))."""
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
    G = expand(rl(b) - rl(b-1))
    return expand(im(G))

def get_Qb(b):
    Ib = ImP_coeff_b(b)
    P = Poly(Ib, m)
    # divide out forced roots m = 0,...,floor((b-1)/2)
    for r in range(0, (b-1)//2 + 1):
        q, rem = sp.div(P, Poly(m - r, m), m)
        assert rem == 0, (b, r, rem)
        P = q
    return P  # this is Q_b (times possible content/sign)

for b in range(1, 21):
    Q = get_Qb(b)
    fac = sp.factor_list(Q.as_expr())
    print(f"b={b:2d} degQ={Q.degree()}  factored: {sp.factor(Q.as_expr())}")
