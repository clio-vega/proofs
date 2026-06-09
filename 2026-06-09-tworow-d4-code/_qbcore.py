"""Reusable core: build primitive integer Q_b for the two-row d=4 law.

I_b(m) = Im G_{(2m-b,b)} = [u^b]((1-u)(1+su+u^2)^m), s=1+i, exact poly in m.
Engine expansion:
  I_b(m) = sum_{j=1, j%4!=0}^{b} C(m,j) Im(s^j) (-1)^{b-j} C(m-j, floor((b-j)/2))
Factor out forced roots prod_{r=0}^{R}(m-r), R=floor((b-1)/2) -> Q_b(m).
Return PRIMITIVE integer Poly (content removed, leading coeff sign positive).
"""
import sympy as sp

m = sp.symbols('m')

def Im_s_pow(j):
    r = j % 4
    if r == 0:
        return 0
    return {1:1, 2:2, 3:2}[r] * (-4)**(j//4)

def _binom_poly(arg, k):
    if k < 0:
        return sp.Integer(0)
    num = sp.Integer(1)
    for i in range(k):
        num *= (arg - i)
    return num / sp.factorial(k)

def I_b(b):
    total = sp.Integer(0)
    for j in range(1, b+1):
        c = Im_s_pow(j)
        if c == 0:
            continue
        eps = (-1)**(b-j)
        k = (b-j)//2
        total += c*eps*_binom_poly(m, j)*_binom_poly(m-j, k)
    return sp.expand(total)

def Qb_primitive(b):
    """Primitive integer Poly Q_b(m), forced roots removed."""
    Ib = I_b(b)
    R = (b-1)//2
    forced = sp.prod([(m-r) for r in range(0, R+1)])
    Qb, rem = sp.div(sp.Poly(Ib, m), sp.Poly(forced, m))
    assert rem == 0, (b, "nonzero remainder")
    Qb = sp.Poly(Qb, m, domain='QQ')
    # clear denominators -> integer, then remove content
    coeffs = Qb.all_coeffs()
    dens = [sp.Rational(c).q for c in coeffs]
    L = sp.ilcm(*dens) if dens else 1
    Qint = sp.Poly([sp.Integer(sp.Rational(c)*L) for c in coeffs], m, domain='ZZ')
    Qprim = Qint.primitive()[1]
    if Qprim.LC() < 0:
        Qprim = -Qprim
    return Qprim

def load_cache(path="results/Qb_cache.pkl"):
    """Return dict b -> sympy Poly Q_b, from the pickled coeff cache."""
    import pickle
    raw = pickle.load(open(path, "rb"))
    return {b: sp.Poly([sp.Integer(c) for c in coeffs], m, domain='ZZ')
            for b, coeffs in raw.items()}

if __name__ == "__main__":
    import time
    for b in [2,3,6,7,30,31,50,51,70,71]:
        t=time.time()
        Q=Qb_primitive(b)
        print(f"b={b:3d} degQ={Q.degree():3d} build={time.time()-t:.2f}s "
              f"LC={Q.LC()} Q(0)={Q.eval(0)}")
