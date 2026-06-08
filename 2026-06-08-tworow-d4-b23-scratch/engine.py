import sympy as sp
m,u = sp.symbols('m u')
s = 1 + sp.I

def g(b, mval):
    P = sp.expand((1 + s*u + u**2)**mval)
    return sp.expand(P).coeff(u,b)

# (a) three-term recurrence:  b*g_b = s(m-b+1)g_{b-1} + (2m-b+2)g_{b-2}
print("=== (a) recurrence check b*g_b = s(m-b+1)g_{b-1}+(2m-b+2)g_{b-2} ===")
ok=True
for mval in range(3,9):
    for b in range(2, 2*mval+1):
        lhs = b*g(b,mval)
        rhs = s*(mval-b+1)*g(b-1,mval) + (2*mval-b+2)*g(b-2,mval)
        if sp.simplify(lhs-rhs)!=0:
            ok=False; print("FAIL",b,mval,sp.simplify(lhs-rhs))
print("recurrence holds:", ok)

# (b) closed form g_b(m) = sum_l falling(m,b-l)/((b-2l)! l!) s^{b-2l}
print("\n=== (b) closed form check ===")
def falling(x,k):
    r=sp.Integer(1)
    for i in range(k): r*= (x-i)
    return r
ok=True
for mval in range(3,9):
    for b in range(0, 2*mval+1):
        cf = sum(falling(mval,b-l)/(sp.factorial(b-2*l)*sp.factorial(l))*s**(b-2*l)
                 for l in range(0,b//2+1) if b-2*l>=0)
        if sp.simplify(cf-g(b,mval))!=0:
            ok=False;print("FAIL",b,mval)
print("closed form holds:", ok)

# (c) Theorem A: I_b((2b-1)/2) = eps_b 2^{floor(b/2)} C(2b,b)/4^b
print("\n=== (c) Theorem A sign/magnitude ===")
from qb import I_b
mm=sp.symbols('m')
for b in range(2,20):
    val=sp.nsimplify(I_b(b).subs(mm, sp.Rational(2*b-1,2)))
    pred_mag = sp.Integer(2)**(b//2)*sp.binomial(2*b,b)/sp.Integer(4)**b
    if val==0:
        print(f"b={b:2d}(mod4={b%4}): 0  (pred mag {pred_mag}, 4|b={b%4==0})")
    else:
        sign = sp.sign(val)
        print(f"b={b:2d}(mod4={b%4}): val={val}, |val|/predmag={sp.nsimplify(abs(val)/pred_mag)}, sign={sign}, "
              f"Im(s^b)sign={sp.sign(sp.im(s**b)) if sp.im(s**b)!=0 else 0}")
