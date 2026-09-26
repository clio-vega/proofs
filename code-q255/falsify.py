"""Q255/Q256: is N(s^c_{lam/mu}) Lorentzian?  Run the falsifier FIRST.

Setup.  h = N(s^c) = sum_alpha (K^c_alpha / alpha!) x^alpha, degree d = |lam/mu|.
LEMMA (recorded, checked below by `hessian_is_raw`):
    the Hessian of d^beta h has entries M_ij = K^c_{beta+e_i+e_j}   -- RAW coeffs,
because d^{beta+e_i+e_j}( (K/gamma!) x^gamma ) = K when gamma = beta+e_i+e_j.
Hence the 2x2-minor consequence of (L3) is, with alpha = beta+e_i+e_j,
    (RLC)   K_alpha^2 >= K_{alpha+e_i-e_j} * K_{alpha-e_i+e_j}
on the RAW coefficients.  The same inequality on the NORMALISED coefficients is a
DIFFERENT and strictly stronger-looking statement; both are reported separately.

All arithmetic is exact (integers / Fraction).
"""
import sys, itertools
from collections import defaultdict
from fractions import Fraction
from math import factorial
sys.path.insert(0, '/home/clio/projects/proofs/code-q254')
sys.path.insert(0, '/home/clio/projects/proofs/code-2026-09-19')
import cyl as C
from winding import weights_counted, winds
from lorentzian import is_M_convex, compositions, multi_factorial


def cyl_shapes(n, m):
    out = []
    for rest in itertools.combinations(range(1, n), m - 1):
        x = (0,) + rest
        if C.is_shape(x, n, m):
            out.append(x)
    return out


def mfact(a):
    r = 1
    for v in a:
        r *= factorial(v)
    return r


# ---------- exact (L3): at most one positive eigenvalue, by sign of char poly ---
def charpoly_int(M):
    """char poly of an integer symmetric matrix by Faddeev-LeVerrier, exact.
    returns ascending [a_0,...,a_n] with det(xI-M) = sum a_k x^k, a_n = 1.
    Recursion (correct indexing):  M_1 = I ;  for k: AM = A M_k,
    c_k = -tr(AM)/k, M_{k+1} = AM + c_k I ;  p(x)=x^n+c_1 x^{n-1}+...+c_n."""
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    Mk = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    cs = []
    for k in range(1, n + 1):
        AM = [[sum(A[i][t] * Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        ck = Fraction(-sum(AM[i][i] for i in range(n)), k)
        cs.append(ck)
        Mk = [[AM[i][j] + (ck if i == j else 0) for j in range(n)] for i in range(n)]
    desc = [Fraction(1)] + cs          # x^n, x^{n-1}, ..., x^0
    return desc[::-1]


def n_positive_eigs_exact(M):
    """# of strictly positive eigenvalues of an integer symmetric matrix, exactly.
    Uses: for a real-rooted poly p(x)=det(xI-M)=sum a_k x^k, the number of positive
    roots equals the number of sign changes in (a_n, a_{n-1}, ..., a_0) restricted to
    nonzero entries (Descartes is EXACT for real-rooted polynomials)."""
    a = charpoly_int(M)                    # ascending
    seq = [c for c in reversed(a) if c != 0]
    return sum(1 for i in range(len(seq) - 1) if (seq[i] > 0) != (seq[i + 1] > 0))


def L3_exact(K, ell, d, want_witness=False):
    """K: dict alpha->int raw coefficient.  Check (L3) exactly via M_ij = K_{b+ei+ej}."""
    if d < 2:
        return True, None
    for beta in compositions(d - 2, ell):
        M = [[0] * ell for _ in range(ell)]
        for i in range(ell):
            for j in range(ell):
                g = list(beta); g[i] += 1; g[j] += 1
                M[i][j] = K.get(tuple(g), 0)
        k = n_positive_eigs_exact(M)
        if k > 1:
            return False, (beta, k, M) if want_witness else (beta, k)
    return True, None


def rlc(K, ell, raw=True):
    """root log-concavity.  raw=True: on K_alpha.  raw=False: on K_alpha/alpha!.
    returns list of violations."""
    bad = []
    for alpha in K:
        for i in range(ell):
            for j in range(ell):
                if i == j:
                    continue
                p = list(alpha); p[i] += 1; p[j] -= 1
                q = list(alpha); q[i] -= 1; q[j] += 1
                if p[j] < 0 or q[i] < 0:
                    continue
                p, q = tuple(p), tuple(q)
                if raw:
                    lhs = K[alpha] ** 2
                    rhs = K.get(p, 0) * K.get(q, 0)
                else:
                    lhs = Fraction(K[alpha], mfact(alpha)) ** 2
                    rhs = (Fraction(K.get(p, 0), mfact(p))
                           * Fraction(K.get(q, 0), mfact(q)))
                if lhs < rhs:
                    bad.append((alpha, i, j, lhs, rhs))
    return bad


def hessian_is_raw():
    """control on the LEMMA above: build N(h) symbolically for a random-ish K and
    confirm d^{beta+ei+ej} N(h) = K_{beta+ei+ej}."""
    import sympy as sp
    ell, d = 3, 4
    xs = sp.symbols('x0 x1 x2')
    K = {}
    for t, alpha in enumerate(compositions(d, ell)):
        K[alpha] = t * t % 7 + 1
    h = sum(sp.Rational(K[a], mfact(a)) * sp.prod([xs[i]**a[i] for i in range(ell)])
            for a in K)
    ok = True
    for beta in compositions(d - 2, ell):
        for i in range(ell):
            for j in range(ell):
                g = list(beta); g[i] += 1; g[j] += 1; g = tuple(g)
                e = h
                for v in range(ell):
                    e = sp.diff(e, xs[v], g[v])
                if sp.simplify(e - K[g]) != 0:
                    ok = False
    return ok
