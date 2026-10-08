"""
2026-10-08 c2 PROVE.  The integral lattice spanned by the power sums.

Engines, deliberately mechanism-disjoint:
  A  M[rho][nu] = #{f:[l(rho)]->[l(nu)] : block sums = nu}      (counting)
  B  M[rho][nu] = [m_nu] p_rho via honest polynomial expansion  (monomials)
  C  M[rho][nu] = z_rho * [p_rho] h_nu                          (Hall pairing in the p-basis)
Engine C never mentions m; engine B never mentions z or h.  Agreement of B and C
is the content of the Q395 gate.
"""
from fractions import Fraction
from itertools import product as iproduct
from functools import lru_cache
from collections import defaultdict

# ---------- partitions ----------
@lru_cache(maxsize=None)
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: return ((),)
    out = []
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            out.append((k,) + rest)
    return tuple(out)

def mults(rho):
    m = defaultdict(int)
    for a in rho: m[a] += 1
    return m

def d_of(rho):
    from math import factorial
    r = 1
    for v in mults(rho).values(): r *= factorial(v)
    return r

def z_of(rho):
    from math import factorial
    r = 1
    for i, v in mults(rho).items(): r *= (i ** v) * factorial(v)
    return r

# ---------- set partitions of [k] ----------
@lru_cache(maxsize=None)
def set_partitions(k):
    if k == 0: return ((),)
    out = []
    for sp in set_partitions(k - 1):
        for i in range(len(sp)):
            out.append(sp[:i] + (sp[i] + (k - 1,),) + sp[i+1:])
        out.append(sp + (((k - 1),),))
    return tuple(out)

def coarsen(lam, pi):
    """lam^pi : block sums, sorted decreasing."""
    return tuple(sorted((sum(lam[i] for i in B) for B in pi), reverse=True))

def mu_hat0(pi):
    """mobius(hat0, pi) in the partition lattice = prod_B (-1)^{|B|-1} (|B|-1)!"""
    from math import factorial
    r = 1
    for B in pi:
        r *= (-1) ** (len(B) - 1) * factorial(len(B) - 1)
    return r

# ---------- Engine A : counting functions ----------
def M_entry_A(rho, nu):
    l, k = len(rho), len(nu)
    cnt = 0
    for f in iproduct(range(k), repeat=l):
        sums = [0] * k
        for i, j in enumerate(f): sums[j] += rho[i]
        if tuple(sums) == tuple(nu): cnt += 1
    return cnt

# ---------- Engine B : honest polynomial expansion in n variables ----------
def p_expand_monomials(rho, nvars):
    """dict: exponent-tuple -> coeff, for p_rho in nvars variables."""
    cur = {tuple([0]*nvars): 1}
    for part in rho:
        nxt = defaultdict(int)
        for exp, c in cur.items():
            for v in range(nvars):
                e = list(exp); e[v] += part
                nxt[tuple(e)] += c
        cur = dict(nxt)
    return cur

def M_row_B(rho, n):
    """row of M from the monomial expansion: coefficient of m_nu = coeff of the
    sorted-exponent representative."""
    nvars = n                      # degree n, at most n parts, so n variables suffice
    poly = p_expand_monomials(rho, nvars)
    row = {}
    for exp, c in poly.items():
        key = tuple(sorted((e for e in exp if e > 0), reverse=True))
        if key in row:
            assert row[key] == c, ("m-symmetry violated", rho, key, row[key], c)
        else:
            row[key] = c
    return row

# ---------- Engine C : Hall pairing via the p-expansion of h ----------
@lru_cache(maxsize=None)
def h_in_p(m):
    """h_m = sum_{sigma |- m} p_sigma / z_sigma  -> dict sigma -> Fraction"""
    return {s: Fraction(1, z_of(s)) for s in partitions(m)}

def hnu_in_p(nu):
    cur = {(): Fraction(1)}
    for part in nu:
        nxt = defaultdict(Fraction)
        for s, c in cur.items():
            for s2, c2 in h_in_p(part).items():
                key = tuple(sorted(s + s2, reverse=True))
                nxt[key] += c * c2
        cur = dict(nxt)
    return cur

def M_col_C(nu):
    """<p_rho, h_nu> = z_rho * [p_rho] h_nu"""
    exp = hnu_in_p(nu)
    return {rho: Fraction(z_of(rho)) * c for rho, c in exp.items()}
