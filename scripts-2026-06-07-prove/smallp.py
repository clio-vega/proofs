import sympy as sp
from sympy import symbols, expand, Poly, I, im
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

def roots_mod(q, p):
    return [r for r in range(p) if q.eval(r) % p == 0]

print("b  b%4  b%8 | roots mod p (p=2,3,5,7) ; '.'=no root (certifies no integer root)")
for b in range(5, 50):
    q = monic_part(b)
    info = []
    nocert = []
    for p in [2,3,5,7,11,13]:
        rs = roots_mod(q, p)
        info.append(f"p{p}:{len(rs)}")
        if len(rs)==0:
            nocert.append(p)
    print(f"{b:2d}  {b%4}   {b%8}  | {' '.join(info)}   no-root-primes={nocert}")
