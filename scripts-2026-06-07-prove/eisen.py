import sympy as sp
from sympy import symbols, expand, Poly, I, im, factorint, gcd
from functools import reduce
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
    Ib = ImP_coeff_b(b)
    P = Poly(Ib, m, domain='QQ')
    for r in range(0, (b-1)//2 + 1):
        P = sp.div(P, Poly(m - r, m, domain='QQ'), m)[0]
    facs = sp.factor_list(P.as_expr())
    big = max((Poly(base, m) for base, mult in facs[1]), key=lambda p: p.degree())
    return Poly(big.as_expr(), m, domain='ZZ')

print("b deg  gcd(non-leading coeffs) factored   | Eisenstein prime candidates (v=1 in c0, divides g)")
for b in range(5, 36):
    q = monic_part(b)
    coeffs = q.all_coeffs()   # leading first
    d = q.degree()
    nonlead = coeffs[1:]      # a_{n-1},...,a_0
    g = reduce(gcd, [abs(int(c)) for c in nonlead])
    c0 = abs(int(coeffs[-1]))
    fg = factorint(g) if g>1 else {}
    # Eisenstein prime: p | g and p^2 does not divide c0
    eis = [p for p in fg if c0 % (p*p) != 0]
    print(f"{b:2d} {d:2d}  g={g} fac={fg}   eisenstein_primes={eis}")
