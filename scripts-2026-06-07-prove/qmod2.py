import sympy as sp
from sympy import symbols, expand, Poly, I, im, GF
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

for b in range(5, 41):
    q = monic_part(b)
    q2 = Poly(q.as_expr(), m, modulus=2)
    fl = q2.factor_list()
    facs = []
    for base, mult in fl[1]:
        facs.append(f"({base.as_expr()})^{mult}" if mult>1 else f"({base.as_expr()})")
    print(f"b={b:2d} b%4={b%4} deg={q.degree()}  q_b mod2 = {' * '.join(facs)}")
