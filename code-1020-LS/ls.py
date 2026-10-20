"""Truncated Laurent-series linear algebra over QQ((t)).

Used to read off ord_t of structure constants at a boundary divisor
{w_m = eps} by substituting w_m = eps + t with all other variables
specialised to random rationals.  Exact (Fraction arithmetic); the only
approximation is the truncation order K, and every routine asserts that it
never needed a coefficient beyond it.
"""
from fractions import Fraction as Fr

K = 40          # default precision: coefficients t^0 .. t^(K-1) after factoring out t^ord


class Trunc(Exception):
    pass


class LS:
    """t^o * (c[0] + c[1] t + ...), c[0] != 0 unless the series is zero."""
    __slots__ = ('o', 'c')

    def __init__(self, o, c):
        # normalise: strip leading zeros
        i = 0
        while i < len(c) and c[i] == 0:
            i += 1
        if i == len(c):
            self.o, self.c = None, []          # zero (to precision)
        else:
            self.o, self.c = o + i, list(c[i:])

    @staticmethod
    def zero():
        return LS(None, [])

    @staticmethod
    def const(x):
        return LS(0, [Fr(x)] + [Fr(0)] * (K - 1))

    def is_zero(self):
        return self.o is None

    def prec(self):
        """number of known coefficients"""
        return len(self.c)

    def __add__(self, other):
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        o = min(self.o, other.o)
        n = min(self.o + len(self.c), other.o + len(other.c)) - o
        if n <= 0:
            raise Trunc("addition: no overlapping known coefficients")
        out = [Fr(0)] * n
        for i, x in enumerate(self.c):
            j = self.o - o + i
            if j < n:
                out[j] += x
        for i, x in enumerate(other.c):
            j = other.o - o + i
            if j < n:
                out[j] += x
        return LS(o, out)

    def __neg__(self):
        return LS(self.o, [-x for x in self.c]) if not self.is_zero() else self

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if self.is_zero() or other.is_zero():
            return LS.zero()
        n = min(len(self.c), len(other.c))
        out = [Fr(0)] * n
        for i, x in enumerate(self.c):
            if i >= n:
                break
            for j, y in enumerate(other.c):
                if i + j >= n:
                    break
                out[i + j] += x * y
        return LS(self.o + other.o, out)

    def inv(self):
        if self.is_zero():
            raise ZeroDivisionError("inverting a series that is zero to precision")
        n = len(self.c)
        a = self.c
        b = [Fr(0)] * n
        b[0] = 1 / a[0]
        for i in range(1, n):
            s = Fr(0)
            for j in range(1, i + 1):
                s += a[j] * b[i - j]
            b[i] = -s / a[0]
        return LS(-self.o, b)

    def __truediv__(self, other):
        return self * other.inv()


def from_poly_coeffs(num, den, prec=None):
    """num, den: lists of Fractions, ascending powers of t.  Return LS of num/den."""
    p = prec or K
    def mk(co):
        return LS(0, (list(co) + [Fr(0)] * p)[:p])
    n, d = mk(num), mk(den)
    if n.is_zero():
        return LS.zero()
    return n / d


def solve(A, b, n):
    """Gauss elimination over QQ((t)) with minimal-valuation pivoting.
    A: n x n list of lists of LS; b: list of LS.  Returns (x, ord_det)."""
    A = [row[:] for row in A]
    b = b[:]
    ordet = 0
    perm_sign_irrelevant = True
    for col in range(n):
        # pivot: smallest valuation
        piv, pv = None, None
        for r in range(col, n):
            e = A[r][col]
            if not e.is_zero() and (pv is None or e.o < pv):
                piv, pv = r, e.o
        if piv is None:
            return None, None          # singular to precision
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            b[col], b[piv] = b[piv], b[col]
        p = A[col][col]
        ordet += p.o
        for r in range(col + 1, n):
            e = A[r][col]
            if e.is_zero():
                continue
            fac = e / p
            for cc in range(col, n):
                A[r][cc] = A[r][cc] - fac * A[col][cc]
            b[r] = b[r] - fac * b[col]
    # back substitution
    x = [LS.zero()] * n
    for r in range(n - 1, -1, -1):
        s = b[r]
        for cc in range(r + 1, n):
            s = s - A[r][cc] * x[cc]
        x[r] = s / A[r][r]
    return x, ordet
