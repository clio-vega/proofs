import sympy as sp
m,u=sp.symbols('m u')
s=1+sp.I

def g_direct(b,mval):
    return sp.expand((1+s*u+u**2)**mval).coeff(u,b)

print("=== verify g_b(m) = C(m,b) s^b 2F1(-b/2,-(b-1)/2; m-b+1; -2i) ===")
ok=True
for b in range(1,9):
    for mval in [5,6,7,9,11]:
        if mval<b: continue
        hyp = sp.binomial(mval,b)*s**b*sp.hyper((-sp.Rational(b,2),-sp.Rational(b-1,2)),(mval-b+1,),-2*sp.I)
        hv = sp.simplify(sp.hyperexpand(hyp) - g_direct(b,mval))
        if hv!=0:
            ok=False; print("FAIL",b,mval, hv)
print("hypergeometric form holds:", ok)

print("\n=== Theorem A as 2F1 at c=1/2 ===")
# g_b((2b-1)/2) = C((2b-1)/2,b) s^b 2F1(-b/2,-(b-1)/2; 1/2; -2i)
for b in range(2,8):
    hyp = sp.binomial(sp.Rational(2*b-1,2),b)*s**b*sp.hyper((-sp.Rational(b,2),-sp.Rational(b-1,2)),(sp.Rational(1,2),),-2*sp.I)
    print(f"b={b}: g_b(half) hyperexpand = {sp.simplify(sp.hyperexpand(hyp))}")
