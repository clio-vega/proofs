"""JOB 2 — exact index map  I_b(m) <-> central trinomial coefficients A002426.

KEY ALGEBRAIC IDENTITY (s = 1+i):
    1 + s u + u^2  =  (1 + u + u^2) + i u  =  W + i u,    W = 1+u+u^2.
Hence
    (1 + s u + u^2)^m = sum_k C(m,k) i^k u^k W^{m-k},
and the imaginary part collects ODD k:
    Im (1+su+u^2)^m = sum_{k odd} (-1)^((k-1)/2) C(m,k) u^k W^{m-k}.
Applying [u^b] (1-u)(.) :
  (TRI)  I_b(m) = sum_{k odd, k<=b} (-1)^((k-1)/2) C(m,k)
                       [ tau(m-k, b-k) - tau(m-k, b-1-k) ],
  tau(n,j) = [u^j] (1+u+u^2)^n  = sum_l C(n,l) C(n-l, j-2l)   (triangle A027907),
  CENTRAL trinomial  T_n := tau(n,n) = A002426 = 1,1,3,7,19,51,141,393,...

So the central trinomial diagonal enters exactly when m-k = b-k, i.e. m = b.
This module:
  (1) verifies W+iu factorisation,
  (2) verifies (TRI) EXACTLY as polynomials in ZZ[m] for b <= 50,
  (3) prints the central-trinomial sequence + its P-recurrence,
  (4) evaluates I_b(b) along the central diagonal and pins the closed shape.
"""
import sympy as sp
from _qbcore import I_b, m

u = sp.symbols('u')
s = 1 + sp.I
W = 1 + u + u**2

def tau_poly(nexpr, j):
    """tau(n,j) = [u^j](1+u+u^2)^n as a polynomial in the symbolic exponent nexpr."""
    if j < 0:
        return sp.Integer(0)
    tot = sp.Integer(0)
    for l in range(0, j//2 + 1):
        a = j - 2*l           # number of u^1 factors
        # C(n,l) * C(n-l, a)
        c1 = sp.Integer(1)
        for i in range(l):
            c1 *= (nexpr - i)
        c1 /= sp.factorial(l)
        c2 = sp.Integer(1)
        for i in range(a):
            c2 *= (nexpr - l - i)
        c2 /= sp.factorial(a)
        tot += c1 * c2
    return sp.expand(tot)

def I_b_trinomial(b):
    """RHS of (TRI) as polynomial in m."""
    tot = sp.Integer(0)
    for k in range(1, b+1, 2):           # odd k
        sign = (-1)**((k-1)//2)
        Cmk = sp.Integer(1)
        for i in range(k):
            Cmk *= (m - i)
        Cmk /= sp.factorial(k)
        tot += sign * Cmk * (tau_poly(m-k, b-k) - tau_poly(m-k, b-1-k))
    return sp.expand(tot)

def central_trinomial(N):
    """T_n = [u^n](1+u+u^2)^n for n=0..N via (n+1)T_{n+1}=(2n+1)T_n+3n T_{n-1}."""
    T = [1, 1]
    for n in range(1, N):
        T.append((( (2*n+1)*T[n] + 3*n*T[n-1] ) // (n+1)))
    return T[:N+1]

if __name__ == "__main__":
    print("=== (1) factorisation 1+su+u^2 = W + i u ===")
    print("  diff =", sp.expand((1 + s*u + u**2) - (W + sp.I*u)))

    print("\n=== (2) EXACT poly identity (TRI) in ZZ[m], b <= 50 ===")
    allok = True
    for b in range(1, 51):
        lhs = I_b(b)
        rhs = I_b_trinomial(b)
        d = sp.expand(lhs - rhs)
        ok = (d == 0)
        allok &= ok
        if not ok:
            print(f"  b={b}: MISMATCH  diff={d}")
    print(f"  (TRI) holds for all b<=50: {allok}")

    print("\n=== (3) central trinomial A002426 + P-recurrence ===")
    T = central_trinomial(12)
    print("  T_n =", T)
    # verify GF 1/sqrt(1-2x-3x^2)
    x = sp.symbols('x')
    ser = sp.series(1/sp.sqrt(1 - 2*x - 3*x**2), x, 0, 9).removeO()
    print("  GF 1/sqrt(1-2x-3x^2) coeffs:",
          [sp.nsimplify(ser.coeff(x, n)) for n in range(9)])

    print("\n=== (4) central diagonal: I_b(b) and its shape ===")
    # along m=b the leading trinomial term is tau(b-k,b-k)=T_{b-k} (central!)
    for b in range(1, 16):
        val = int(I_b(b).subs(m, b))
        print(f"  b={b:2d}: I_b(b) = {val}")
