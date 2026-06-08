import sympy as sp
from qb import I_b, binom_poly, Im_s_pow
m = sp.symbols('m')

def Qprim_of(b):
    Ib = I_b(b)
    R = (b-1)//2
    forced = sp.prod([(m-r) for r in range(0, R+1)])
    Qb, rem = sp.div(sp.Poly(Ib, m), sp.Poly(forced, m))
    assert rem == 0
    coeffs = sp.Poly(Qb, m).all_coeffs()
    dens = [sp.fraction(c)[1] for c in coeffs]
    L = sp.lcm_list(dens) if dens else 1
    Qint = [sp.nsimplify(c*L) for c in coeffs]
    cont = sp.gcd(Qint)
    return sp.Poly([c//cont for c in Qint], m)

print("=== explicit Q_b (primitive) ===")
for b in [2,3,6,7,10,11]:
    Q = Qprim_of(b)
    cc=Q.all_coeffs()
    print(f"b={b}: Q_b = {Q.as_expr()}")
    print(f"      const={cc[-1]} factored={sp.factorint(int(cc[-1])) if cc[-1]!=0 else 0}, lead={cc[0]}")

print("\n=== half-integer (2b-1)/2 eval ===")
for b in range(2,31):
    if b%4 not in (2,3): continue
    Ib = I_b(b)
    val = sp.nsimplify(Ib.subs(m, sp.Rational(2*b-1,2)))
    print(f"b={b:2d}: I_b((2b-1)/2) = {val}")

print("\n=== signs at m=b..b+25 ===")
for b in [6,7,10,11,14,15,18,19]:
    Ib = I_b(b)
    s=''.join('+' if int(Ib.subs(m,mv))>0 else ('-' if int(Ib.subs(m,mv))<0 else '0') for mv in range(b,b+26))
    print(f"b={b:2d}: {s}")
