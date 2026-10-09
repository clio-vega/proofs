"""Theorem D' : the diagonal family.  All claims measured exactly where possible."""
import sympy as sp, mpmath as mp, sys
t = sp.Symbol('t')
mp.mp.dps = 40

def D(b): return sp.Poly(t**b - t**(b-1) + 1, t)
def E(b): return sp.Poly(t**b - t**(b-1) + 2, t)

def circle_roots(poly, tol=mp.mpf('1e-25')):
    """roots of modulus 1, as mpmath complex numbers (high precision)."""
    cs = [mp.mpf(int(c)) for c in poly.all_coeffs()]
    rts = mp.polyroots(cs, maxsteps=400, extraprec=400)
    return [r for r in rts if abs(abs(r)-1) < tol], rts

def mahler_from_roots(poly):
    cs = [mp.mpf(int(c)) for c in poly.all_coeffs()]
    rts = mp.polyroots(cs, maxsteps=400, extraprec=400)
    lead = abs(cs[0])
    return lead*mp.fprod([max(mp.mpf(1), abs(r)) for r in rts]), rts

def mahler_from_jensen(poly):
    """M(f) = exp( int_0^1 log|f(e^{2 pi i th})| dth )  -- Jensen's formula, a DIFFERENT mechanism."""
    cs = [int(c) for c in poly.all_coeffs()]
    def f(th):
        z = mp.e**(2j*mp.pi*th)
        v = mp.mpf(0)
        acc = mp.mpc(0)
        for c in cs: acc = acc*z + c
        return mp.log(abs(acc))
    return mp.e**mp.quad(f, [0,1])

def nu(poly):
    """nu(f) = prod over nonzero roots of max(|a|,1/|a|)."""
    cs = [mp.mpf(int(c)) for c in poly.all_coeffs()]
    rts = mp.polyroots(cs, maxsteps=400, extraprec=400)
    out = mp.mpf(1)
    for r in rts:
        a = abs(r)
        if a > mp.mpf('1e-30'): out *= max(a, 1/a)
    return out

Phi6 = sp.Poly(t**2 - t + 1, t)

if __name__ == '__main__':
    BMAX = int(sys.argv[1]) if len(sys.argv)>1 else 40
    print("b | deg | #roots on |z|=1 | those roots | Phi6 | simple | M(roots) | M(Jensen) | nu | nu-M^2")
    bad_circle=0; bad_phi=0; bad_simple=0; bad_nu=0; bad_cross=0; checked=0
    for b in range(1, BMAX+1):
        p = D(b); checked += 1
        onc, rts = circle_roots(p)
        q, r = sp.div(p, Phi6, t)
        phi6_div = (sp.expand(r.as_expr()) == 0)
        simple = sp.gcd(p, p.diff(t)).degree() == 0
        Mr, _ = mahler_from_roots(p); Mj = mahler_from_jensen(p); nv = nu(p)
        # CLAIM A: roots on |z|=1 are exactly the two primitive 6th roots of unity, and
        #          they occur iff b = 2 mod 6
        want = 2 if (b % 6 == 2) else 0
        if len(onc) != want: bad_circle += 1
        if onc:
            for r0 in onc:
                if abs(r0**6 - 1) > mp.mpf('1e-20') or abs(r0**3-1) < mp.mpf('1e-20'):
                    bad_circle += 1
        # CLAIM B: Phi6 | D_b  iff  b = 2 mod 6   (exact division, not numerics)
        if phi6_div != (b % 6 == 2): bad_phi += 1
        # CLAIM C: all roots simple
        if not simple: bad_simple += 1
        # CLAIM D: nu(D_b) > 1 for b >= 3 ; = 1 for b = 1,2
        if b >= 3 and not (nv > 1 + mp.mpf('1e-6')): bad_nu += 1
        if b <= 2 and abs(nv - 1) > mp.mpf('1e-20'): bad_nu += 1
        # CROSS-CHECK: two mechanisms for M agree, and nu = M^2 (valid since |D_b(0)|=1)
        if abs(Mr - Mj) > mp.mpf('1e-12') or abs(nv - Mr**2) > mp.mpf('1e-12'): bad_cross += 1
        if b <= 12 or b % 6 == 2 or b == BMAX:
            print("%2d | %3d | %d | %s | %s | %s | %s | %s | %s | %.2e" % (
                b, p.degree(), len(onc),
                ",".join(mp.nstr(r,6) for r in onc) if onc else "-",
                "yes" if phi6_div else "no", "yes" if simple else "NO",
                mp.nstr(Mr,10), mp.nstr(Mj,10), mp.nstr(nv,10), float(abs(nv-Mr**2))))
    print()
    print("D_b, b=1..%d : %d values checked" % (BMAX, checked))
    print("  CLAIM A  unit-circle roots are exactly {e^{+-i pi/3}} and only when b=2 mod 6 : %d failures" % bad_circle)
    print("  CLAIM B  Phi6 | D_b  <=>  b = 2 mod 6  (exact polynomial division)           : %d failures" % bad_phi)
    print("  CLAIM C  every root of D_b is simple                                         : %d failures" % bad_simple)
    print("  CLAIM D  nu(D_b) > 1 for b>=3, nu(D_b) = 1 for b in {1,2}                    : %d failures" % bad_nu)
    print("  CROSS    M(roots) = M(Jensen) and nu = M^2                                   : %d failures" % bad_cross)
