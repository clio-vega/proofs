"""
TASK 3: Identification sweep (FAST version).
Match I_b(m) and primitive Q_b(m) up to rational scalar and affine reparam
m -> alpha*m + beta against classical OP families.

Speed: precompute each family poly ONCE per param set, get its rational coeff
vector via numpy-free exact arithmetic, substitute x->alpha*m+beta numerically-
exactly via sympy Poly.transform, then compare *monic coefficient vectors*.
"""
import sys
import sympy as sp
from sympy import symbols, Rational, Poly, expand, rf, binomial, Integer
from identify_sweep import deflate, primitive_int_poly, I_b, m as M

x = symbols('x')

def polyexpand(e):
    """Force binomial(var,k) etc into explicit polynomials."""
    return expand(sp.expand_func(e))

def monic_vec(poly_expr, var):
    p = Poly(polyexpand(poly_expr), var)
    lc = p.LC()
    return tuple(c/lc for c in p.all_coeffs()), p.degree()

def poch(a, k):
    return rf(a, k)

def meixner(n, beta, c):
    tot = Integer(0)
    for k in range(n+1):
        tot += poch(-n,k)*poch(-x,k)/poch(beta,k)*(1-Rational(1,1)/c)**k/sp.factorial(k)
    return expand(tot)

def krawtchouk_q(n, N, q):
    tot = Integer(0)
    for j in range(n+1):
        tot += (q-1)**(n-j)*(-1)**j * binomial(x,j)*binomial(N-x, n-j)
    return expand(tot)

def hahn(n, alpha, beta, N):
    tot = Integer(0)
    for k in range(n+1):
        denom = poch(alpha+1,k)*poch(-N,k)*sp.factorial(k)
        if denom == 0: return None
        tot += poch(-n,k)*poch(n+alpha+beta+1,k)*poch(-x,k)/denom
    return expand(tot)

def dual_hahn(n, gamma, delta, N):
    tot = Integer(0)
    for k in range(n+1):
        denom = poch(gamma+1,k)*poch(-N,k)*sp.factorial(k)
        if denom == 0: return None
        tot += poch(-n,k)*poch(-x,k)*poch(x+gamma+delta+1,k)/denom
    return expand(tot)

def disc_cheb(n, N):
    tot = Integer(0)
    for k in range(n+1):
        tot += (-1)**k * binomial(n,k)*binomial(n+k,k)*binomial(x,k)
    return expand(tot)

ALPHAS = (Integer(1), Integer(-1), Integer(2), Rational(1,2))
SHIFTS = range(-3,4)

def fam_monic_vecs(fam_expr, Tdeg):
    """Return dict mapping monic-vec -> (alpha,beta) over reparams, if degree matches."""
    fam_p = polyexpand(fam_expr)
    p = Poly(fam_p, x)
    if p.degree() != Tdeg:
        return {}
    out = {}
    for alpha in ALPHAS:
        for beta in SHIFTS:
            sub = Poly(expand(fam_p.subs(x, alpha*M + Integer(beta))), M)
            if sub.degree() != Tdeg:
                continue
            lc = sub.LC()
            vec = tuple(cc/lc for cc in sub.all_coeffs())
            out.setdefault(vec, (alpha, beta))
    return out

def scan(label, fam_func, params, target_vec, Tdeg):
    matches = []
    for prm in params:
        fam = fam_func(Tdeg, *prm)
        if fam is None:
            continue
        try:
            vecs = fam_monic_vecs(fam, Tdeg)
        except Exception:
            continue
        if target_vec in vecs:
            a,bb = vecs[target_vec]
            matches.append((prm, a, bb))
    return matches

def main():
    bs = [5,6,7,8,9,10]
    betas = [Rational(1,2),1,Rational(3,2),2,3]
    cs = [Rational(1,2),-1,2,Rational(-1,2),Rational(1,3),3,-2]
    Ns = list(range(2,18))
    qs = [Rational(1,2),2,3,-1,Rational(3,2),4]

    for b in bs:
        Ib, Qraw, forced = deflate(b)
        coeffs, Qprim = primitive_int_poly(Qraw)
        for tname, target in [("I_b", Ib), ("Q_b", Qprim)]:
            tvec, Tdeg = monic_vec(target, M)
            print(f"\n--- b={b} target={tname} deg={Tdeg} ---"); sys.stdout.flush()
            any_found=False

            mm = scan("Meixner", meixner, [(bt,c) for bt in betas for c in cs], tvec, Tdeg)
            if mm: any_found=True; print(f"  MEIXNER (beta,c),alpha,shift: {mm[:8]}")

            mm = scan("Kraw", lambda n,N,q: krawtchouk_q(n,N,q),
                      [(N,q) for N in Ns for q in qs], tvec, Tdeg)
            if mm: any_found=True; print(f"  KRAWTCHOUK-q (N,q),alpha,shift: {mm[:8]}")

            mm = scan("Hahn", hahn,
                      [(al,bt,N) for al in betas for bt in betas for N in Ns if N>=Tdeg],
                      tvec, Tdeg)
            if mm: any_found=True; print(f"  HAHN (alpha,beta,N),alpha,shift: {mm[:8]}")

            mm = scan("dHahn", dual_hahn,
                      [(g,d,N) for g in betas for d in betas for N in Ns if N>=Tdeg],
                      tvec, Tdeg)
            if mm: any_found=True; print(f"  DUAL-HAHN (gamma,delta,N),alpha,shift: {mm[:8]}")

            mm = scan("dCheb", lambda n,N: disc_cheb(n,N),
                      [(N,) for N in Ns if N>=Tdeg], tvec, Tdeg)
            if mm: any_found=True; print(f"  DISCRETE-CHEBYSHEV (N,),alpha,shift: {mm[:8]}")

            if not any_found:
                print("  no match in scanned families/params")
            sys.stdout.flush()

if __name__=="__main__":
    main()
