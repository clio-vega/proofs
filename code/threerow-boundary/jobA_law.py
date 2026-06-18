"""
jobA_law.py -- consolidate the exact content law, indexed by (c,i) not just depth k.

Computes the fixed-divisor content of master N_i^{(c)} on each slice (grid known
exact from jobA_verify_content), organizes by i and by parity of i, and tests:
  H1 (FLOOR, the certificate PROVE relies on): content >= 2*floor(k/2) always;
     even-slice content >= k always.
  H2 (i-parity law): on the ODD slice, ODD i  => content == 2*floor(k/2) (floor);
                                      EVEN i => content > floor (surplus).
  Also reports content - v2(c!k!) and the residual beyond floor, to expose the
  c-dependence that defeats any (k, c mod 4) closed form.
"""
import sympy as sp
from math import factorial
from fast_alt import Malt

P, B = sp.symbols('P B')


def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def Ni_val(c, i, Pv, Bv, even):
    k = c - i
    bv = 2 * Bv if even else 2 * Bv + 1
    av = 2 * Pv + bv + c
    m = Pv + bv + c
    j = bv + i
    if bv < c + 1 or j > m:
        return None
    M = Malt(av, bv, c, j, m)
    den = bv - c + 1
    for t in range(2, i + 1):
        den *= (bv + t)
    num = M * factorial(c) * factorial(k)
    if num % den:
        return None
    return num // den


def fit(c, i, even, pdeg, bdeg):
    monos = [(p, q) for p in range(pdeg + 1) for q in range(bdeg + 1)]
    B0 = ((c + 1) // 2) + 2
    pts, rhs = [], []
    for Pv in range(pdeg + 1):
        for Bv in range(B0, B0 + bdeg + 1):
            v = Ni_val(c, i, Pv, Bv, even)
            if v is None:
                return None
            pts.append((Pv, Bv)); rhs.append(sp.Integer(v))
    A = sp.Matrix([[sp.Integer(Pv)**p * sp.Integer(Bv)**q for (p, q) in monos]
                   for (Pv, Bv) in pts])
    sol = A.LUsolve(sp.Matrix(rhs))
    return sp.expand(sum(sol[t]*P**p*B**q for t, (p, q) in enumerate(monos)))


def content(expr, even, c, grid=96):
    poly = sp.Poly(sp.expand(expr), P, B)
    D = 1
    for cf in poly.coeffs():
        D = sp.ilcm(D, sp.Rational(cf).q)
    vD = v2int(D)
    terms = [(int(m[0]), int(m[1]), int(sp.Rational(cf)*D)) for m, cf in poly.terms()]
    qmax = max(t[1] for t in terms)
    g = 10**9
    B0 = ((c+1)//2)+1
    for Bv in range(B0, grid):
        Bp = [Bv**q for q in range(qmax+1)]
        for Pv in range(grid):
            val = sum(cc*(Pv**pp)*Bp[qq] for (pp, qq, cc) in terms)
            vv = v2int(val) - vD
            if vv < g:
                g = vv
    return g


if __name__ == "__main__":
    data = {}   # (c,i,even) -> content
    for c in range(4, 13):
        for i in range(2, c+1):
            k = c-i
            if k > 6:
                continue
            pdeg, bdeg = k//2+3, k+4
            for even in (True, False):
                e = fit(c, i, even, pdeg, bdeg)
                data[(c, i, even)] = content(e, even, c)
        print("c=%d done" % c)

    print("\n" + "="*90)
    print("CONTENT by (c,i) : even / odd  ;  k=c-i ; floor=2*floor(k/2)")
    print("="*90)
    h1_floor = h1_evenk = True
    for c in range(4, 13):
        for i in range(2, c+1):
            k = c-i
            if k > 6:
                continue
            ge, go = data[(c, i, True)], data[(c, i, False)]
            floor = 2*(k//2)
            if ge < floor or go < floor:
                h1_floor = False
            if ge < k:
                h1_evenk = False
            print("  c=%-2d i=%-2d k=%-2d (i %s) | even=%d odd=%d | floor=%d  vk=%d  %s%s"
                  % (c, i, k, "odd " if i % 2 else "even", ge, go, floor, k,
                     "EVEN<k!" if ge < k else "",
                     " ODDfloor" if (i % 2 == 1 and go == floor) else ""))

    print("\nH1 FLOOR (content >= 2floor(k/2) both slices):", h1_floor)
    print("H1 even-slice >= k:", h1_evenk)

    # H2: odd-i odd-slice == floor ?
    print("\n" + "="*90)
    print("H2: ODD slice, by parity of i")
    print("="*90)
    for ip in (1, 0):
        print(" i %s:" % ("ODD" if ip else "EVEN"))
        for k in range(0, 7):
            vals = sorted(set(data[(c, c-k, False)] - 2*(k//2)
                              for c in range(4, 13)
                              if (c, c-k, False) in data and (c-k) % 2 == ip and c-k >= 2))
            print("   k=%d : odd-slice surplus over floor = %s" % (k, vals))
    print("\n   (surplus 0 == content sits exactly at the floor)")

    # i-parity on even slice surplus over k
    print("\nEVEN slice surplus over k (content - k), by parity of i:")
    for ip in (1, 0):
        print(" i %s:" % ("ODD" if ip else "EVEN"))
        for k in range(0, 7):
            vals = sorted(set(data[(c, c-k, True)] - k
                              for c in range(4, 13)
                              if (c, c-k, True) in data and (c-k) % 2 == ip and c-k >= 2))
            print("   k=%d : even-slice (content - k) = %s" % (k, vals))
