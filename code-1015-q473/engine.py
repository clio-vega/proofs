"""
Q473 engine: lim_{c->infty} G_kappa(x_1..x_L | Y) of Garbali-de Gier-Mead-Wheeler,
arXiv:2312.14348, Theorem `G nu thm` (line 3386 of Half_line_ASEP.tex).

Pure SymPy. No SageMath (verified absent 2026-10-06).

Conventions, measured off the source:
  h(x)        = a c (1-x^2) / ((a-x)(c-x))                     [eq h-intro, line 279]
  Q(xi,xj)    = (1-h(xi))(1-h(xj)) - h(xi)h(xj)/(ac) * (1-q)xi xj/(1-q xi xj)
                                                               [eq: Q Pfaffian kernel, l.540]
  Z_m         = prod_{i<j} (1-xi xj)/(xi-xj) * Pf( (xi-xj)/(1-xi xj) Q(xi,xj) )   m even
  Z_{2l-1}(x_1..x_{2l-1}) := Z_{2l}(x_1,..,x_{2l-1},1)         [remark after thm:ZPfaff]
  lim_{c->oo} Z_L = prod_i (1-h^oo(x_i)) = prod_i x_i(1-a x_i)/(x_i-a)   [line ~620]
  phi_m(x)    = (1-q) x y_m/(1-q x y_m) * prod_{j=1}^{m-1} (1-x y_j)/(1-q x y_j)
                 (the paper prints `\prod_{j-1}^{\nu_i-1}` in two later displays; in
                  `G nu sum expression thm eq` itself it correctly prints `\prod_{j=1}`.)
"""
import sympy as sp
from itertools import combinations, permutations

q, a, c = sp.symbols('q a c')

def ysym(j):
    return sp.Symbol('y%d' % j)

def xsym(i):
    return sp.Symbol('x%d' % i)

# ---------- generic c ----------
def h(x):
    return a*c*(1-x**2)/((a-x)*(c-x))

def Qker(xi, xj):
    return (1-h(xi))*(1-h(xj)) - h(xi)*h(xj)/(a*c) * (1-q)*xi*xj/(1-q*xi*xj)

def pfaffian(M):
    """Pfaffian of an even-dimensional antisymmetric sympy Matrix, by recursive expansion
    along the first row:  Pf(M) = sum_{j>1} (-1)^j M[0,j] Pf(M with rows/cols 0,j removed)."""
    n = M.shape[0]
    if n == 0:
        return sp.Integer(1)
    assert n % 2 == 0
    tot = sp.Integer(0)
    for j in range(1, n):
        idx = [k for k in range(1, n) if k != j]
        sub = M[idx, idx]
        tot += (-1)**(j+1) * M[0, j] * pfaffian(sub)
    return tot

def Z_generic(xs):
    """Triangular partition function Z_m at generic boundary parameters a, c."""
    xs = list(xs)
    m = len(xs)
    if m == 0:
        return sp.Integer(1)
    if m % 2 == 1:
        return Z_generic(xs + [sp.Integer(1)])
    M = sp.zeros(m, m)
    for i in range(m):
        for j in range(m):
            if i != j:
                M[i, j] = (xs[i]-xs[j])/(1-xs[i]*xs[j]) * Qker(xs[i], xs[j])
    pref = sp.Integer(1)
    for i in range(m):
        for j in range(i+1, m):
            pref *= (1-xs[i]*xs[j])/(xs[i]-xs[j])
    return sp.together(sp.cancel(sp.together(pref*pfaffian(M))))

# ---------- c -> infinity ----------
def h_inf(x):
    return a*(1-x**2)/(a-x)

def Z_inf(xs):
    r = sp.Integer(1)
    for x in xs:
        r *= x*(1-a*x)/(x-a)
    return r

def phi(m, x):
    """phi_m(x); m >= 1."""
    r = (1-q)*x*ysym(m)/(1-q*x*ysym(m))
    for j in range(1, m):
        r *= (1-x*ysym(j))/(1-q*x*ysym(j))
    return r

def G(kappa, L, Zfun=Z_inf, hfun=h_inf, sign_plant=False):
    """G_kappa(x_1..x_L) via `G nu sum expression thm eq`.
    kappa: tuple, strictly DECREASING (kappa_1 > ... > kappa_d >= 1).
    Zfun/hfun let the caller swap in the generic-c versions.
    sign_plant: planted error control -- flips a sign in the S_n product."""
    kappa = tuple(kappa)
    d = len(kappa)
    assert L >= d, "G nu thm requires L >= n"
    assert all(kappa[i] > kappa[i+1] for i in range(d-1)), "kappa must be strictly decreasing"
    xs = [xsym(i) for i in range(1, L+1)]
    total = sp.Integer(0)
    for K in combinations(range(L), d):
        Kb = [i for i in range(L) if i not in K]
        term = Zfun([xs[i] for i in Kb])
        for i in K:
            term *= hfun(xs[i])
        for i in K:
            for j in Kb:
                term *= (xs[j]-q*xs[i])/(xs[j]-xs[i]) * (1-xs[i]*xs[j])/(1-q*xs[i]*xs[j])
        for ii in range(d):
            for jj in range(ii+1, d):
                term *= (1-xs[K[ii]]*xs[K[jj]])/(1-q*xs[K[ii]]*xs[K[jj]])
        inner = sp.Integer(0)
        for sigma in permutations(range(d)):
            t = sp.Integer(1)
            for ii in range(d):
                for jj in range(ii+1, d):
                    num = xs[K[sigma[jj]]] - q*xs[K[sigma[ii]]]
                    if sign_plant:
                        num = xs[K[sigma[jj]]] + q*xs[K[sigma[ii]]]
                    t *= num/(xs[K[sigma[jj]]] - xs[K[sigma[ii]]])
            for ii in range(d):
                t *= phi(kappa[ii], xs[K[sigma[ii]]])
            inner += t
        total += term*inner
    return total

def subsets_upto(maxel, maxcard):
    """all strictly-decreasing tuples from {1..maxel} of cardinality <= maxcard."""
    out = []
    for d in range(0, maxcard+1):
        for s in combinations(range(1, maxel+1), d):
            out.append(tuple(reversed(s)))
    return out

def dominates(kp, k):
    """kp >= k in the order of Lemma (triangularity): |kp|>=|k| and kp_i >= k_i for i<=|k|."""
    if len(kp) < len(k):
        return False
    return all(kp[i] >= k[i] for i in range(len(k)))
