"""
Two-row d=4 law, b = 2,3 mod 4.  Study Q_b(m).

I_b(m) = Im G_{(2m-b,b)} = [u^{b-1}]( (1-u) H_m ),  exact integer poly in m.
Exact expansion (engine):
  I_b(m) = sum_{j=1, j !=0 mod4}^{b}  C(m,j) Im(s^j) (-1)^{b-j} C(m-j, floor((b-j)/2)),
  s=1+i, Im(s^j) = (-4)^{floor(j/4)} * {1,2,2}[j mod4 in {1,2,3}].
Factor out forced roots prod_{r=0}^{R}(m-r), R=floor((b-1)/2)  -> Q_b(m).
"""
import sympy as sp

m = sp.symbols('m')

def Im_s_pow(j):
    r = j % 4
    if r == 0:
        return 0
    base = {1:1, 2:2, 3:2}[r]
    return base * (-4)**(j//4)

def binom_poly(arg, k):
    """C(arg, k) as explicit polynomial in m: prod_{i=0}^{k-1}(arg-i)/k!."""
    if k < 0:
        return sp.Integer(0)
    num = sp.Integer(1)
    for i in range(k):
        num *= (arg - i)
    return sp.expand(num / sp.factorial(k))

def I_b(b):
    """Return I_b(m) as sympy poly in m via exact engine expansion."""
    total = sp.Integer(0)
    for j in range(1, b+1):
        c = Im_s_pow(j)
        if c == 0:
            continue
        eps = (-1)**(b-j)
        k = (b-j)//2
        binom1 = binom_poly(m, j)
        binom2 = binom_poly(m-j, k)
        total += c*eps*binom1*binom2
    return sp.expand(total)

def I_b_direct(b, mval):
    """Direct check: Im G via [u^b]((1-u)(1+su+u^2)^m) using gaussian ints."""
    u = sp.symbols('u')
    s = 1 + sp.I
    P = sp.expand((1 - u)*(1 + s*u + u**2)**mval)
    coeff = P.coeff(u, b)
    return sp.im(coeff)

# sanity check the engine vs direct
print("=== sanity: engine vs direct ===")
for b in [2,3,6,7,10,11]:
    Ib = I_b(b)
    for mval in range(b, b+4):
        eng = int(Ib.subs(m, mval))
        dr = int(I_b_direct(b, mval))
        assert eng == dr, (b, mval, eng, dr)
    print(f"b={b}: OK")

print()
print("=== Q_b structure for b = 2,3 mod 4 ===")
for b in range(2, 31):
    if b % 4 not in (2,3):
        continue
    Ib = I_b(b)
    R = (b-1)//2
    forced = sp.prod([(m-r) for r in range(0, R+1)])
    Qb, rem = sp.div(sp.Poly(Ib, m), sp.Poly(forced, m))
    assert rem == 0, (b, "nonzero remainder!")
    Qb = sp.Poly(Qb, m)
    # primitive integer version
    coeffs = Qb.all_coeffs()
    dens = [sp.fraction(sp.nsimplify(c))[1] for c in coeffs]
    L = sp.lcm_list(dens) if dens else 1
    Qint = sp.Poly([sp.nsimplify(c*L) for c in coeffs], m)
    cont = sp.gcd(list(Qint.all_coeffs()))
    Qprim = sp.Poly([c//cont for c in Qint.all_coeffs()], m)
    fact = sp.factor_list(Qprim.as_expr())
    disc = sp.discriminant(Qprim.as_expr(), m)
    # real roots >= b ?
    realroots = sp.Poly(Qprim, m).real_roots()
    rr = [sp.nsimplify(r) for r in realroots]
    rr_num = sorted(float(r) for r in realroots)
    big = [round(x,3) for x in rr_num if x >= b-0.5]
    irr = "irred" if len(fact[1])==1 and fact[1][0][1]==1 else "REDUCIBLE"
    print(f"b={b:2d} (mod4={b%4}) degQ={Qprim.degree():2d} {irr:10s} "
          f"realroots>=b: {big}")
    if irr != "irred":
        print("        factors:", fact)
