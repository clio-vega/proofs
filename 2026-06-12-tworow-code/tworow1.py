import sympy as sp

x, u = sp.symbols('x u')

# Conjectured Prop-1 (two-row):
#   Sum_j C(m,j) M_j x^j = ([u^a] - [u^{a+1}]) (u^2 + (2+x)u + 1)^m
# Derived closed form claim:
#   M_j = C(2(m-j), b-j) - C(2(m-j), b-j-1)   (using b = min row, b<=m)
# i.e. M_j = f^{(a-j,b-j)}.
#
# We test against the DEFINITION M_j = <s_lambda, h1^{2(m-j)} e2^j>.
# But computing that needs symmetric functions; here first just test
# the internal consistency: does extracting [x^j] of the boxed GF give the binomial form?

def poly_coeffs(m, a):
    g = (u**2 + (2+x)*u + 1)**m
    g = sp.expand(g)
    g = sp.Poly(g, u)
    # [u^a] - [u^{a+1}]
    ca = g.coeff_monomial(u**a) if a>=0 else 0
    ca1 = g.coeff_monomial(u**(a+1))
    expr = sp.expand(ca - ca1)
    return expr  # polynomial in x = Sum_j C(m,j) M_j x^j

def Mj_from_gf(m,a,b,j):
    expr = poly_coeffs(m,a)
    px = sp.Poly(expr, x)
    cj = px.coeff_monomial(x**j)   # = C(m,j) M_j
    return sp.Integer(cj)

from sympy import binomial as C

ok=True
for m in range(1,8):
    for b in range(0, m+1):
        a = 2*m - b
        for j in range(0, b+1):
            lhs = Mj_from_gf(m,a,b,j)            # = C(m,j) M_j
            # closed form using b-coordinates:
            Mj = C(2*(m-j), b-j) - C(2*(m-j), b-j-1)
            rhs = C(m,j)*Mj
            if sp.simplify(lhs-rhs)!=0:
                ok=False
                print("MISMATCH", m,a,b,j, lhs, rhs)
print("closed-form vs GF consistent:", ok)
