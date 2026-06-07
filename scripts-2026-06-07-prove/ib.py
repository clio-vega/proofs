import sympy as sp

m, u = sp.symbols('m u', real=True)
I = sp.I

def Ib_poly(b, Mmax=None):
    """Return I_b(m) as a sympy polynomial in m.
    I_b(m) = Im [u^{b}]((1-u)(1+su+u^2)^m), s=1+i.
    We compute via [u^l]P = sum_c multinom(m; l-2c,c,m-l+c) s^{l-2c}, l=b and l=b-1?
    Careful: I_b(m) = [u^b]((1-u) Im P) where Im P. Actually from writeup:
    G_{(2m-b,b)} = [u^b]((1-u)P), I_b = Im of that = [u^b]((1-u) Im P)...
    but Im P only; (1-u) real. So I_b(m) = [u^b](Im P) - [u^{b-1}](Im P).
    """
    # build symbolic m-polynomial for [u^l] Im P
    def coeff_ImP(l):
        # sum over c>=0 with l-2c>=0 and m-(l-c)>=0 (handle as polynomial in m, factorials -> falling)
        total = sp.Integer(0)
        c = 0
        while l - 2*c >= 0:
            a = l - 2*c
            # multinomial m!/(a! c! (m-a-c)!) = C(m,a+c)*C(a+c,c)*... as poly in m:
            # = (1/(a! c!)) * m*(m-1)*...*(m-(a+c)+1)
            ac = a + c
            fall = sp.Integer(1)
            for j in range(ac):
                fall *= (m - j)
            multinom = fall / (sp.factorial(a)*sp.factorial(c))
            s_pow = (1+I)**a
            total += multinom * s_pow
            c += 1
        total = sp.expand(total)
        # take imaginary part: coefficients are Gaussian integers times m-poly
        return sp.im(total)  # this treats m as real symbol; sp.im on expr with I
    Cl_b = coeff_ImP(b)
    Cl_bm1 = coeff_ImP(b-1)
    res = sp.expand(Cl_b - Cl_bm1)
    return sp.Poly(res, m)

for b in range(1, 13):
    P = Ib_poly(b)
    print(f"b={b}: deg={P.degree()}")
    print("   I_b(m) =", sp.factor(P.as_expr()))
    # rational roots
    rr = sp.roots(P, m)
    print("   roots:", {sp.nsimplify(k): v for k,v in rr.items()})
