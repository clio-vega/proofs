import sympy as sp
m,u=sp.symbols('m u')
s=1+sp.I
def G_half(b):
    mm=sp.Rational(2*b-1,2)
    P=sp.expand((1-u)*(1+s*u+u**2)**mm)  # fractional power -> sympy series
    # use closed form for g_b - g_{b-1} via coefficient of u^b in (1-u)A^m
    # compute via the verified engine closed form instead:
    def falling(x,k):
        r=sp.Integer(1)
        for i in range(k): r*=(x-i)
        return r
    def g(b):
        return sum(falling(mm,b-l)/(sp.factorial(b-2*l)*sp.factorial(l))*s**(b-2*l)
                   for l in range(0,b//2+1) if b-2*l>=0)
    return sp.expand(g(b)-g(b-1))
print("=== verify  G_b((2b-1)/2) = (-1)^b C(2b,b)/4^b (1-i)^b ===")
ok=True
for b in range(1,16):
    lhs=sp.simplify(G_half(b))
    rhs=sp.simplify((-1)**b*sp.binomial(2*b,b)/sp.Integer(4)**b*(1-sp.I)**b)
    d=sp.simplify(lhs-rhs)
    print(f"b={b:2d}: G_b(half)={lhs}   match={d==0}")
    if d!=0: ok=False
print("CLOSED FORM HOLDS:", ok)
print("\n=== therefore I_b(half)=Im G = (-1)^(b+1) C(2b,b)/4^b Im(s^b) ===")
for b in range(2,12):
    lhs=sp.im(sp.simplify(G_half(b)))
    rhs=(-1)**(b+1)*sp.binomial(2*b,b)/sp.Integer(4)**b*sp.im(s**b)
    print(f"b={b:2d}: Im={lhs}, formula={sp.nsimplify(rhs)}, match={sp.simplify(lhs-rhs)==0}, |.|=2^{{floor(b/2)}}C/4^b? {sp.simplify(abs(lhs)-sp.Integer(2)**(b//2)*sp.binomial(2*b,b)/sp.Integer(4)**b)==0 if b%4!=0 else 'N/A(4|b ->0)'}")
