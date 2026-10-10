"""Q405 core: symmetric-function engine for Warnaar Eq_generalise vs HKKO Thm 3.3.

Conventions, all read verbatim from source (2026-10-10):

HKKO 2301.13117 Thm 3.3 (thm:C), eq (C-):
    det_{1<=i,j<=h}( F_{-i+j,2h+w}(x) - F_{i+j,2h+w}(x) )
      = sum_{la in Par(2h,w)} c^-_{2h,w}(la) * s_{la[2h,w]'}(x)
with
    f_r(x)       = sum_{i in Z} e_i(x) e_{i+r}(x)          (f_r = f_{-r})
    F_{r,N}(x)   = sum_{k in Z} f_{r+Nk}(x)
    Fbar_{r,N}(x)= sum_{k in Z} (-1)^k f_{r+Nk}(x)
    c^-_{2h,w}(la) = +1 if la_{2i-1}=la_{2i} for all 1<=i<=h
                     -1 if la_1-la_{2h}=w and la_{2i}=la_{2i+1} for 1<=i<=h-1
                      0 otherwise
    Par(2h,w) = { la : l(la) <= 2h, la_1 - la_{2h} <= w }

HKKO Thm 2.x (eq:CJT), cylindric Jacobi-Trudi:
    s_{la[m,w]'}(x) = sum_{k in Z^m, sum k_i = 0} det_{1<=i,j<=m}( e_{la_i-i+j+(m+w)k_i}(x) )

Warnaar 2511.17034 (Eq_generalise), with M := k+l+1, N := 2M, h := k (det size),
w := 2l+2 so that N = 2h + w:
    G_t = sum_{y in Z^h} det_{1<=i,j<=h}( t^{M y_i^2 - j y_i}
              ( edot_{n-i+j-N y_i}(x) - edot_{n+i+j-N y_i}(x) ) )
    edot_r(x) = e_r(x_1,1/x_1,...,x_n,1/x_n)
    Warnaar: "Eq_generalise for t=1 is (x_1...x_n)^{-k} det(F_{i-j,N} - F_{i+j,N})"
"""
import itertools, functools
from sympy import symbols, Poly, expand, simplify, Rational, together, cancel, factor
import sympy as sp


def elem(xs):
    """e_0..e_n in the variables xs (list)."""
    n = len(xs)
    e = [sp.Integer(0)]*(n+1)
    e[0] = sp.Integer(1)
    for x in xs:
        for m in range(n, 0, -1):
            e[m] = sp.expand(e[m] + x*e[m-1])
    return e


class Engine:
    def __init__(self, n):
        self.n = n
        self.x = sp.symbols(f'x1:{n+1}', positive=True)
        self.e = elem(list(self.x))                      # e_m(x), 0<=m<=n
        self.prodx = sp.prod(self.x)
        # edot_r = e_r(x^{\pm}), 0<=r<=2n
        self.edot = elem([v for x in self.x for v in (x, 1/x)])

    def E(self, m):
        return self.e[m] if 0 <= m <= self.n else sp.Integer(0)

    def Edot(self, r):
        return self.edot[r] if 0 <= r <= 2*self.n else sp.Integer(0)

    @functools.lru_cache(maxsize=None)
    def f(self, r):
        """f_r = sum_i e_i e_{i+r}; f_r = f_{-r}."""
        r = abs(r)
        if r > self.n:
            return sp.Integer(0)
        return sp.expand(sum(self.E(i)*self.E(i+r) for i in range(0, self.n+1)))

    @functools.lru_cache(maxsize=None)
    def F(self, r, N):
        tot = sp.Integer(0)
        k = 0
        while True:
            added = False
            for rr in ({r + N*k, r - N*k} if k else {r}):
                if abs(rr) <= self.n:
                    tot += self.f(rr); added = True
            if k > 0 and not added and abs(r) - N*k < -self.n and r + N*k > self.n:
                break
            k += 1
            if N*k > 2*self.n + abs(r) + 4:
                break
        return sp.expand(tot)
