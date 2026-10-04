"""Exact-arithmetic Lorentzian / M-convexity toolkit for (Q), 2026-10-04 c3.

All arithmetic is exact integer arithmetic.  No floats anywhere in a signature test.

Conventions
-----------
A homogeneous polynomial  f = sum_alpha c_alpha v^alpha   (alpha in Z_{>=0}^n, |alpha|=d)
is stored as a dict  {alpha (tuple) : c_alpha (int)}  with c_alpha > 0 on the support.
`N(f)` is never formed explicitly: by the raw-Hessian lemma the Hessian of
d^beta N(f) is the integer matrix [c_{beta+e_i+e_j}], so every Lorentzian test runs
on the RAW coefficients c.
"""
from itertools import product as iproduct
from fractions import Fraction

# ---------------------------------------------------------------- M-convexity

def is_Mconvex(S):
    """Symmetric exchange axiom, brute force.  S a set/list of integer tuples."""
    S = set(S)
    if not S:
        return True
    n = len(next(iter(S)))
    sums = {sum(a) for a in S}
    if len(sums) != 1:
        return False                      # M-convex => constant coordinate sum
    for a in S:
        for b in S:
            for i in range(n):
                if a[i] > b[i]:
                    ok = False
                    for j in range(n):
                        if a[j] < b[j]:
                            a2 = list(a); a2[i] -= 1; a2[j] += 1
                            b2 = list(b); b2[i] += 1; b2[j] -= 1
                            if tuple(a2) in S and tuple(b2) in S:
                                ok = True; break
                    if not ok:
                        return False
    return True

# ------------------------------------------------- positive-eigenvalue counts

def _char_poly_sym3(M):
    """(e1,e2,e3) of a symmetric 3x3 integer matrix."""
    e1 = M[0][0] + M[1][1] + M[2][2]
    e2 = (M[0][0]*M[1][1] - M[0][1]**2
          + M[0][0]*M[2][2] - M[0][2]**2
          + M[1][1]*M[2][2] - M[1][2]**2)
    e3 = (M[0][0]*(M[1][1]*M[2][2] - M[1][2]**2)
          - M[0][1]*(M[0][1]*M[2][2] - M[1][2]*M[0][2])
          + M[0][2]*(M[0][1]*M[1][2] - M[1][1]*M[0][2]))
    return e1, e2, e3

def at_most_one_pos_3x3_detcrit(M):
    """l3-det-reduction (lean-verified): for symmetric 3x3 with NONNEGATIVE entries,
    at most one positive eigenvalue  <=>  e_2 <= 0 and det >= 0."""
    e1, e2, e3 = _char_poly_sym3(M)
    return e2 <= 0 and e3 >= 0

def count_pos_eigen_exact(M):
    """Number of strictly positive eigenvalues of a symmetric integer matrix,
    WITH MULTIPLICITY, computed exactly.

    DEFECT FOUND AND FIXED 2026-10-04 c3.  The first version stripped the root at 0
    and returned `q.count_roots(0, oo)`.  sympy's count_roots counts DISTINCT real
    roots in the interval, so the function undercounted repeated eigenvalues: it
    returned 1 for M = 4*I_2, whose eigenvalue 4 has multiplicity 2.  Since this
    function was the INDEPENDENT validator of the lean-verified criterion
    l3-det-reduction, the defect sat upstream of that validation.  It was caught by a
    single apparent violation of Lemma minor, which was a false alarm about the
    lemma and a true alarm about the instrument.

    `Poly.real_roots()` returns real roots WITH multiplicity, and a real symmetric
    matrix has only real eigenvalues, so the list has length n."""
    import sympy as sp
    lam = sp.Symbol('lam')
    p = sp.Poly(sp.Matrix(M).charpoly(lam).as_expr(), lam)
    roots = p.real_roots()
    assert len(roots) == len(M), "symmetric matrix must have n real eigenvalues"
    return sum(1 for r in roots if r > 0)

def at_most_one_pos(M, use_sympy=False):
    if len(M) == 3 and not use_sympy and all(x >= 0 for r in M for x in r):
        return at_most_one_pos_3x3_detcrit(M)
    return count_pos_eigen_exact(M) <= 1

# ------------------------------------------------- the Lorentzian test on N(f)

def hessian(c, beta, n):
    """[c_{beta+e_i+e_j}]_{i,j}  -- the raw-Hessian lemma."""
    M = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            g = list(beta); g[i] += 1; g[j] += 1
            M[i][j] = c.get(tuple(g), 0)
    return M

def betas(d, n):
    """all beta in Z_{>=0}^n with |beta| = d-2"""
    if d < 2:
        return []
    def rec(k, rem):
        if k == 1:
            yield (rem,); return
        for v in range(rem+1):
            for t in rec(k-1, rem-v):
                yield (v,)+t
    return list(rec(n, d-2))

def is_N_lorentzian(c, check_support=True, use_sympy=False, return_witness=False):
    """Is N(f) Lorentzian, for f with raw coefficient dict c?
       Criterion (Branden-Huh \\label{SecondDefinition} l.637 + \\label{chars} l.1623):
       support M-convex, and for every |beta|=d-2 the Hessian of d^beta N(f) has at
       most one positive eigenvalue."""
    supp = [a for a, v in c.items() if v != 0]
    if not supp:
        return (True, None) if return_witness else True
    n = len(supp[0])
    d = sum(supp[0])
    if any(sum(a) != d for a in supp):
        return (False, 'inhomogeneous') if return_witness else False
    if any(v < 0 for v in c.values()):
        return (False, 'negative coefficient') if return_witness else False
    if check_support and not is_Mconvex(supp):
        return (False, 'support not M-convex') if return_witness else False
    for b in betas(d, n):
        M = hessian(c, b, n)
        if not at_most_one_pos(M, use_sympy=use_sympy):
            return (False, ('beta', b, M)) if return_witness else False
    return (True, None) if return_witness else True

# ------------------------------------------------------------- the polynomials

def poly_mult(f, g):
    h = {}
    for a, u in f.items():
        for b, v in g.items():
            k = tuple(x+y for x, y in zip(a, b))
            h[k] = h.get(k, 0) + u*v
    return h

def Ptilde(l, h):
    """P~(X,E,W) = sum_{x+e+w=h, l<=x+e<=h} X^x E^e W^w.  All coefficients 1.
       Support = {(x,e,w)>=0 : x+e+w=h, w <= h-l}."""
    c = {}
    for w in range(0, h-l+1):
        for x in range(0, h-w+1):
            c[(x, h-w-x, w)] = 1
    return c

def k_sequence(L, H, D):
    """k(a) = #{(x,e) in Z_{>=0}^{2m} : L_i<=x_i+e_i<=H_i, sum x=a, sum e=D-a},
       by direct enumeration (independent of the Lorentzian machinery)."""
    m = len(L)
    # dp over beads: state = (sum x, sum (x+e))
    cur = {(0, 0): 1}
    for i in range(m):
        nxt = {}
        for (sx, su), v in cur.items():
            for u in range(L[i], H[i]+1):
                for x in range(0, u+1):
                    key = (sx+x, su+u)
                    nxt[key] = nxt.get(key, 0) + v
        cur = nxt
    out = {}
    for (sx, su), v in cur.items():
        if su == D:
            out[sx] = out.get(sx, 0) + v
    return [out.get(a, 0) for a in range(0, D+1)]

def is_logconcave_pf2(seq):
    """PF_2: log-concave coefficients AND no internal zeros."""
    nz = [i for i, v in enumerate(seq) if v != 0]
    if not nz:
        return True
    lo, hi = nz[0], nz[-1]
    if any(seq[i] == 0 for i in range(lo, hi+1)):
        return False
    for a in range(lo+1, hi):
        if seq[a]*seq[a] < seq[a-1]*seq[a+1]:
            return False
    return True
