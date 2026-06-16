import sympy as sp
a,b,P = sp.symbols('a b P', integer=True)

# c=3: a = 2P+b+3
a_sub = 2*P+b+3

# N2 = a(b+3) - (b^2+2b+3)
N2 = a*(b+3) - (b**2+2*b+3)
N2s = sp.expand(N2.subs(a,a_sub))
print("N2 in (P,b):", sp.factor(N2s))
# my claim: N2 = 2[(P+1)(b+3)+b]
mine = 2*((P+1)*(b+3)+b)
print("N2 - mine =", sp.simplify(N2s-mine))

# N1 = a b^2 + 5ab + 6a - b^3 - b^2 - 4b - 6
N1 = a*b**2+5*a*b+6*a - b**3-b**2-4*b-6
N1s = sp.expand(N1.subs(a,a_sub))
print("N1 in (P,b):", sp.factor(N1s))
