import sympy as sp
from qb import I_b
m = sp.symbols('m')

print("=== I_b((2b-1)/2) vs central binomial C(2b,b)/4^b ===")
for b in range(2,24):
    Ib = I_b(b)
    val = sp.nsimplify(Ib.subs(m, sp.Rational(2*b-1,2)))
    cb = sp.binomial(2*b,b)/sp.Integer(4)**b
    ratio = sp.nsimplify(val/cb)
    print(f"b={b:2d}(mod4={b%4}): I_b(half)={str(val):22s} ratio to C(2b,b)/4^b = {ratio} = {sp.nsimplify(ratio)}")

print("\n=== general: I_b at m = half-integer  (2t+1)/2 for several t ===")
# evaluate I_b at m=(2t-1)/2 and compare to central binomials
for b in [6,7,10,11]:
    Ib=I_b(b)
    print(f"b={b}:")
    for t in range(b, b+6):
        val=sp.nsimplify(Ib.subs(m, sp.Rational(2*t-1,2)))
        # odd part / 2-power
        print(f"   m=(2*{t}-1)/2={sp.Rational(2*t-1,2)}: {val}")
