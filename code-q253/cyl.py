"""Cylindric skew Schur functions from the CSSYT definition (AKO Def 2.x / def:cylindricSkewSchur).

Conventions verified against AKO 2311.07382 examples (7-box CSSYT example, ex:mnrule):
  * cylinder C_{x,y}: (0,0) ~ (x,-y) in Z^2, i.e. box (r,c) ~ (r+y, c+x)   [r = drawn row, down]
  * a cylindric diagram is given by row intervals [lo_i, hi_i], i in Z,
    with lo_{i+y} = lo_i + x, hi_{i+y} = hi_i + x, lo/hi weakly increasing in i.
  * CSSYT: entries weakly increase as column index increases (rightward)
           and strictly increase as row index DEcreases (upward in the drawing).
  * loop ribbon = ribbon with x+y boxes; height of a loop ribbon is y,
    height of a non-loop ribbon is (#rows it meets) - 1.
"""
from itertools import product
from functools import lru_cache

class Cyl:
    def __init__(self, x, y, lo, hi):
        assert len(lo) == len(hi) == y
        self.x, self.y = x, y
        self.lo = list(lo); self.hi = list(hi)
        self.n = x + y
        for i in range(y):
            assert self.LO(i) <= self.LO(i+1), ("lo not weakly incr", lo)
            assert self.HI(i) <= self.HI(i+1), ("hi not weakly incr", hi)
            assert self.HI(i) >= self.LO(i)-1
    def LO(self, i):
        q, r = divmod(i, self.y); return self.lo[r] + q*self.x
    def HI(self, i):
        q, r = divmod(i, self.y); return self.hi[r] + q*self.x
    def cells(self):
        "representative cells: rows 0..y-1"
        return [(i,c) for i in range(self.y) for c in range(self.LO(i), self.HI(i)+1)]
    def size(self):
        return len(self.cells())
    def rep(self, i, c):
        "canonical representative of cell (i,c)"
        q, r = divmod(i, self.y); return (r, c - q*self.x)
    def inside(self, i, c):
        return self.LO(i) <= c <= self.HI(i)
    def __repr__(self):
        return f"Cyl(x={self.x},y={self.y},lo={self.lo},hi={self.hi},|D|={self.size()})"

def cssyt_monomials(D, nvars):
    """dict: sorted weight tuple -> count.  Brute force over fillings."""
    cells = D.cells()
    idx = {c:k for k,c in enumerate(cells)}
    N = len(cells)
    # constraint list on representative indices
    weak = []   # (a,b) : T[a] <= T[b]
    strict = [] # (a,b) : T[a] <  T[b]
    for (i,c) in cells:
        if D.inside(i, c+1):
            weak.append((idx[(i,c)], idx[D.rep(i,c+1)]))
        # vertical: (i,c) above (i+1,c) -> T(i,c) > T(i+1,c)
        if D.inside(i+1, c):
            strict.append((idx[D.rep(i+1,c)], idx[(i,c)]))
    from collections import Counter
    out = Counter()
    T = [0]*N
    def rec(k):
        if k == N:
            cnt = Counter(T)
            out[tuple(cnt.get(v,0) for v in range(1,nvars+1))] += 1
            return
        for v in range(1, nvars+1):
            T[k] = v
            ok = True
            for (a,b) in weak:
                if a<=k and b<=k and T[a] > T[b]: ok=False; break
            if ok:
                for (a,b) in strict:
                    if a<=k and b<=k and T[a] >= T[b]: ok=False; break
            if ok: rec(k+1)
        T[k]=0
    rec(0)
    return dict(out)

# ---- symmetric function bookkeeping -------------------------------------
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k): yield (k,)+rest

def zee(mu):
    from collections import Counter
    from math import factorial
    c = Counter(mu); z = 1
    for part, m in c.items(): z *= part**m * factorial(m)
    return z

def mono_coeffs(D, nvars=None):
    """coefficient of m_lambda in s_D, for all partitions lambda of |D| with <= nvars parts."""
    N = D.size()
    if nvars is None: nvars = N
    raw = cssyt_monomials(D, nvars)
    out = {}
    for lam in partitions(N):
        if len(lam) > nvars: continue
        key = tuple(lam) + (0,)*(nvars-len(lam))
        out[lam] = raw.get(key, 0)
    return out

def m_to_p(N, nvars=None):
    """expansion of m_lambda in the p basis is awkward; instead go m -> monomial eval.
    We use:  <s, p_nu> = c_nu  where s = sum c_nu p_nu / z_nu.
    Easiest route: expand s in the *power sum* basis by first writing s in terms of
    monomials and using the transition via h/e is messy.  Instead we use
    the augmented-monomial <-> p duality:  p_nu = sum_lambda R_{nu,lambda} m_lambda
    where R_{nu,lambda} = #{ ways to distribute the parts of nu into blocks with sizes lambda }.
    Then c_nu is obtained by solving the triangular system.
    """
    raise NotImplementedError

def p_mono_coeff(nu, lam):
    """coefficient of the monomial x_1^{lam_1} ... x_k^{lam_k} in p_nu = prod_i (sum_j x_j^{nu_i}).
    = # assignments of the parts of nu to the slots 1..k with prescribed slot sums."""
    k = len(lam)
    from functools import lru_cache
    nu = tuple(nu)
    @lru_cache(maxsize=None)
    def rec(i, rem):
        if i == len(nu):
            return 1 if all(r == 0 for r in rem) else 0
        tot = 0
        for j in range(k):
            if rem[j] >= nu[i]:
                nr = list(rem); nr[j] -= nu[i]
                tot += rec(i+1, tuple(nr))
        return tot
    return rec(0, tuple(lam))

def p_expansion(mcoeffs, N):
    """given {lambda: coeff of m_lambda}, return {nu: c_nu} with s = sum c_nu p_nu / z_nu."""
    from fractions import Fraction
    parts = list(partitions(N))
    A = [[Fraction(p_mono_coeff(nu, lam), zee(nu)) for nu in parts] for lam in parts]
    b = [Fraction(mcoeffs.get(lam, 0)) for lam in parts]
    # gaussian elimination
    n = len(parts)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f*bb for a, bb in zip(M[r], M[col])]
    return {parts[i]: M[i][n] for i in range(n)}

def conjugate(D):
    """transpose (r,c) -> (c,r).  D on C_{x,y} -> D' on C_{y,x}.
    new row c = [B(c), A(c)] with A(c)=max{r: LO(r)<=c}, B(c)=min{r: HI(r)>=c}
    (empty when B(c) = A(c)+1)."""
    x, y = D.x, D.y
    R = range(-6*(x+y)-6, 6*(x+y)+6)
    def A(c): return max(r for r in R if D.LO(r) <= c)
    def B(c): return min(r for r in R if D.HI(r) >= c)
    c0 = min(c for (_, c) in D.cells())
    nlo = [B(c) for c in range(c0, c0+x)]
    nhi = [A(c) for c in range(c0, c0+x)]
    sh = nlo[0]
    return Cyl(y, x, [v-sh for v in nlo], [v-sh for v in nhi])

def eps(nu):
    return (-1)**(sum(nu) - len(nu))

def c_nu(D, nvars=None):
    mc = mono_coeffs(D, nvars=nvars if nvars else D.size())
    return p_expansion(mc, D.size())

def enumerate_cyl(x, y, maxsize=None):
    """all cylindric diagrams on C_{x,y} up to column translation (lo_0 = 0)."""
    out = []
    def los(i, cur):
        if i == y:
            yield tuple(cur); return
        for v in range(cur[-1] if cur else 0, x+1):
            if i == 0 and v != 0: continue
            yield from los(i+1, cur+[v])
    for lo in los(0, []):
        if lo[-1] > lo[0] + x: continue
        def his(i, cur):
            if i == y:
                yield tuple(cur); return
            start = max(lo[i]-1, cur[-1] if cur else lo[i]-1)
            top = (cur[0] + x) if cur else (lo[i]-1 + x + y)
            for v in range(start, top+1):
                yield from his(i+1, cur+[v])
        for hi in his(0, []):
            if hi[-1] > hi[0] + x: continue
            try: D = Cyl(x, y, lo, hi)
            except AssertionError: continue
            if maxsize is not None and D.size() > maxsize: continue
            if D.size() == 0: continue
            out.append(D)
    return out
