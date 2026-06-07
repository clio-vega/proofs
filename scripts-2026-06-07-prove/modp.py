import sympy as sp
from sympy import symbols, expand, Poly, I, im, factorint, GF
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

print("b deg  smallest p with q_b irreducible mod p  | #irred-factors-mod-3 | const-factorization")
for b in range(5, 36):
    q = monic_part(b)
    d = q.degree()
    irr = sp.factor_list(q.as_expr())
    n_irr = len(irr[1])  # over Q: should be 1 if irreducible
    # smallest prime making it irreducible mod p
    pp = None
    for p in sp.primerange(3, 200):
        try:
            fp = Poly(q.as_expr(), m, modulus=p)
            fl = fp.factor_list()
            if len(fl[1]) == 1 and fl[1][0][1] == 1 and fl[1][0][0].degree() == d:
                pp = p
                break
        except Exception:
            continue
    c0 = q.all_coeffs()[-1]
    print(f"{b:2d} {d:2d}  irr/Q={n_irr==1}  p_irr={pp}  const={factorint(abs(int(c0)))}")
