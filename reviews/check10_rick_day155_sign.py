"""Rick Day-155 (commit 88194e3) section 1 boxes, at n=3:
        Psi^+(f)(u) = - Psi(f o phi)(-u),      phi: u_i -> -u_i
   using V(-u) = (-1)^{binom(n,2)} V(u) = -V(u) at n=3.

   Test it directly from HIS OWN definitions Psi(f) = T(fV)/V, Psi^+(f) = T^+(fV)/V,
   T: u^alpha -> prod (u_i)_{alpha_i} (falling), T^+: u^alpha -> prod u_i^{(alpha_i)} (rising).
   Uses a NON-homogeneous f as well, so homogeneity cannot mask a global sign.
"""
import sympy as sp

def falling(y,k):
    e=sp.Integer(1)
    for t in range(k): e*=(y-t)
    return e
def rising(y,k):
    e=sp.Integer(1)
    for t in range(k): e*=(y+t)
    return e

def apply_T(poly, u, ff):
    """Apply T (or T^+) to a polynomial: monomialwise u^alpha -> prod ff(u_i, alpha_i)."""
    p = sp.Poly(sp.expand(poly), *u)
    out = sp.Integer(0)
    for alpha, c in zip(p.monoms(), p.coeffs()):
        term = sp.Integer(c)
        for i, a in enumerate(alpha): term *= ff(u[i], a)
        out += term
    return sp.expand(out)

def V(u):
    n=len(u)
    return sp.expand(sp.prod([u[i]-u[j] for i in range(n) for j in range(i+1,n)]))

def schur(mu, u):
    n=len(u); k=[(mu[j] if j<len(mu) else 0)+n-1-j for j in range(n)]
    return sp.cancel(sp.expand(sp.Matrix(n,n,lambda i,j: u[i]**k[j]).det())/V(u))

def Psi(f, u, ff):
    return sp.cancel(apply_T(sp.expand(f*V(u)), u, ff)/V(u))

for n in (3,4):
    u = sp.symbols(f'u0:{n}')
    phi = {ui: -ui for ui in u}
    tests = {'s_(2,1) (homogeneous)': schur((2,1),u),
             's_(2,1) + s_(1) (NON-homogeneous)': schur((2,1),u)+schur((1,),u)}
    print(f"n={n}  (binom(n,2) = {n*(n-1)//2}, so V(-u) = {(-1)**(n*(n-1)//2):+d} V(u))")
    for name, f in tests.items():
        lhs = sp.expand(Psi(f, u, rising))
        f_phi = sp.expand(f.subs(phi, simultaneous=True))
        rhs_plus = sp.expand(Psi(f_phi, u, falling).subs(phi, simultaneous=True))
        eq_plus  = sp.simplify(lhs - rhs_plus) == 0
        eq_minus = sp.simplify(lhs + rhs_plus) == 0
        print(f"   {name:36s}  Psi^+ == +Psi(f.phi)(-u): {eq_plus};   == -Psi(f.phi)(-u): {eq_minus}")
