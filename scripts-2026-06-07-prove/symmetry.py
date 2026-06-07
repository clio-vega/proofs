import sympy as sp
m = sp.symbols('m', real=True)
I = sp.I

def coeff_ImP(l):
    total = sp.Integer(0); c = 0
    while l - 2*c >= 0:
        a = l - 2*c; ac = a + c
        fall = sp.Integer(1)
        for j in range(ac): fall *= (m - j)
        total += fall/(sp.factorial(a)*sp.factorial(c)) * (1+I)**a
        c += 1
    return sp.im(sp.expand(total))

def Ib(b):
    return sp.expand(coeff_ImP(b) - coeff_ImP(b-1))

# Search for reflection symmetry I_b(m) = eps * I_b(c - m), eps=+-1, c rational with small denom
print("=== Reflection symmetry search: I_b(m) =? eps * I_b(c-m) ===")
for b in range(1, 13):
    P = Ib(b)
    found = []
    for eps in (1, -1):
        # solve for c: P(m) - eps P(c-m) == 0 as identity in m. Try c in halves.
        for c2 in range(-2, 8*b+4):  # c = c2/2
            c = sp.Rational(c2, 2)
            diff = sp.expand(P - eps*P.subs(m, c - m))
            if diff == 0:
                found.append((eps, c))
    print(f"b={b}: deg={sp.Poly(P,m).degree()}  symmetries (eps,c): {found}")
