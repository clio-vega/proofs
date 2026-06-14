"""
Core library: SYT enumeration, Des/maj/s(T), G_lambda(q), fake degrees,
2-core / 2-quotient, and verification utilities.

Self-contained. Uses sympy for exact polynomial arithmetic.
"""
import sympy as sp
from sympy import symbols, Poly, factor_list, prod, factorial, binomial, cyclotomic_poly
from functools import lru_cache
from itertools import product as iproduct

q = symbols('q')


# ---------------------------------------------------------------------------
# Partitions
# ---------------------------------------------------------------------------
def partitions(n):
    """All partitions of n as weakly-decreasing tuples."""
    if n == 0:
        return [()]
    res = []
    def gen(remaining, maxpart, cur):
        if remaining == 0:
            res.append(tuple(cur))
            return
        for p in range(min(remaining, maxpart), 0, -1):
            gen(remaining - p, p, cur + [p])
    gen(n, n, [])
    return res


# ---------------------------------------------------------------------------
# SYT enumeration
# ---------------------------------------------------------------------------
def syt_list(lam):
    """
    Enumerate standard Young tableaux of shape lam.
    Returns list of fillings: filling[r][c] = entry (1..n).
    """
    lam = list(lam)
    n = sum(lam)
    rows = len(lam)
    # cells in row-major; we place 1..n
    cells = [(r, c) for r in range(rows) for c in range(lam[r])]
    result = []

    # represent tableau as list of lists, -1 empty
    tab = [[0] * lam[r] for r in range(rows)]
    for r in range(rows):
        for c in range(lam[r]):
            tab[r][c] = -1

    def addable_cells(filled_count_per_row):
        # cell (r,c) is addable if c == current length filled in row r,
        # and (r==0 or filled in row r-1 > c)
        out = []
        for r in range(rows):
            c = filled_count_per_row[r]
            if c < lam[r]:
                if r == 0 or filled_count_per_row[r - 1] > c:
                    out.append(r)
        return out

    filled = [0] * rows

    def rec(k):
        if k > n:
            # record
            snapshot = [row[:] for row in tab]
            result.append(snapshot)
            return
        for r in addable_cells(filled):
            c = filled[r]
            tab[r][c] = k
            filled[r] += 1
            rec(k + 1)
            filled[r] -= 1
            tab[r][c] = -1

    rec(1)
    return result


def descent_set(tab):
    """Des(T) = {i : i+1 in strictly lower row than i}."""
    # map entry -> row
    row_of = {}
    for r, row in enumerate(tab):
        for v in row:
            row_of[v] = r
    n = len(row_of)
    des = set()
    for i in range(1, n):
        if row_of[i + 1] > row_of[i]:
            des.add(i)
    return des


def maj(tab):
    return sum(descent_set(tab))


def s_stat(tab, n):
    """s(T) = sum_{i in Des} w_i, w_i = 2i-1 if (n-i) odd else 0."""
    des = descent_set(tab)
    total = 0
    for i in des:
        if (n - i) % 2 == 1:
            total += 2 * i - 1
    return total


@lru_cache(maxsize=None)
def G_lambda(lam):
    """G_lambda(q) = sum_T q^{s(T)} as sympy Poly in q."""
    lam = tuple(lam)
    n = sum(lam)
    tabs = syt_list(lam)
    poly = sp.Integer(0)
    for t in tabs:
        poly += q ** s_stat(t, n)
    return sp.expand(poly)


@lru_cache(maxsize=None)
def fake_degree(mu):
    """f^mu(q) = sum_{U in SYT(mu)} q^{maj(U)}."""
    mu = tuple(mu)
    if sum(mu) == 0:
        return sp.Integer(1)
    tabs = syt_list(mu)
    poly = sp.Integer(0)
    for t in tabs:
        poly += q ** maj(t)
    return sp.expand(poly)


def num_syt(lam):
    return len(syt_list(lam))


# ---------------------------------------------------------------------------
# 2-core and 2-quotient via beta-sets / abacus
# ---------------------------------------------------------------------------
def beta_set(lam, r):
    """
    First-column hook lengths / beta-numbers: beta_i = lam_i + (r-1-i) for i=0..r-1,
    using r beads (r >= len(lam)). Returns list of length r (strictly decreasing).
    """
    lam = list(lam) + [0] * (r - len(lam))
    return [lam[i] + (r - 1 - i) for i in range(r)]


def two_core_and_quotient(lam):
    """
    Compute the 2-core and 2-quotient of lam using the abacus (2 runners).
    Returns (core_partition, (q0, q1)) where q0,q1 are partitions.
    Standard James/Kerber convention.
    """
    lam = list(lam)
    # choose number of beads r even, r >= len(lam), large enough
    r = len(lam) + (len(lam) % 2)  # make even
    r = max(r, 2)
    if r % 2 != 0:
        r += 1
    # ensure enough beads
    while r < len(lam):
        r += 2
    beta = beta_set(lam, r)
    beta_set_vals = set(beta)

    # Distribute beads onto 2 runners by residue mod 2.
    # Position p -> runner (p mod 2), row = p // 2.
    runners = {0: [], 1: []}
    for b in beta:
        runners[b % 2].append(b // 2)
    for k in runners:
        runners[k].sort(reverse=True)

    # quotient: for runner k with beads at positions (// 2) values v_0>v_1>..,
    # number of beads m_k, the partition is v_i - (m_k-1-i)
    quotient = []
    for k in (0, 1):
        vals = runners[k]
        m = len(vals)
        part = [vals[i] - (m - 1 - i) for i in range(m)]
        part = [x for x in part if x > 0]
        quotient.append(tuple(part))

    # 2-core: push all beads up on each runner (remove gaps), reconstruct partition
    core_beta = []
    for k in (0, 1):
        m = len(runners[k])
        # pushed-up positions: 0,1,...,m-1 on runner k => abacus positions k + 2*row
        for row in range(m):
            core_beta.append(k + 2 * row)
    core_beta.sort(reverse=True)
    rc = len(core_beta)
    core = [core_beta[i] - (rc - 1 - i) for i in range(rc)]
    core = tuple(x for x in core if x > 0)

    return core, (quotient[0], quotient[1])


def n_of_lambda(lam):
    """n(lam) = sum_i (i-1) lam_i  (1-indexed rows)."""
    return sum(i * part for i, part in enumerate(lam))  # i from 0 => (i-1) with 1-index = i here


def chi_w0_via_domino(lam):
    """
    Predicted chi^lambda(w0) via the 2-quotient domino formula:
      0 unless 2-core is empty (n even) or single cell (n odd).
      else (-1)^{n(lam)} * f^{q0} * f^{q1} * binom(floor(n/2), |q0|).
    """
    core, (q0, q1) = two_core_and_quotient(lam)
    n = sum(lam)
    core_size = sum(core)
    if n % 2 == 0:
        if core_size != 0:
            return 0
    else:
        if core_size != 1:
            return 0
    f0 = num_syt(q0)
    f1 = num_syt(q1)
    return ((-1) ** n_of_lambda(lam)) * f0 * f1 * int(binomial(n // 2, sum(q0)))


# ---------------------------------------------------------------------------
# Cyclotomic factorization
# ---------------------------------------------------------------------------
def cyclotomic_factorization(poly):
    """
    Return dict {d: mult} so poly = leading * prod Phi_d(q)^mult, plus leftover.
    Returns (content, {d:mult}, leftover_poly).
    Strategy: divide out cyclotomics up to degree.
    """
    p = sp.Poly(sp.expand(poly), q)
    if p.is_zero:
        return (0, {}, sp.Integer(0))
    factors = {}
    deg = p.degree()
    # try cyclotomic Phi_d for d up to (deg+1)*2 to be safe
    maxd = 2 * (deg + 2) + 4
    changed = True
    while changed:
        changed = False
        for d in range(1, maxd + 1):
            phi = sp.Poly(cyclotomic_poly(d, q), q)
            if phi.degree() == 0:
                continue
            quo, rem = sp.div(p, phi)
            if rem.is_zero and quo.degree() >= 0:
                factors[d] = factors.get(d, 0) + 1
                p = quo
                changed = True
                if p.degree() == 0:
                    break
    leftover = p.as_expr()
    content = leftover  # whatever remains (should be constant if fully cyclotomic)
    return (content, factors, leftover)


def fmt_cyclo(factors):
    if not factors:
        return "1"
    parts = []
    for d in sorted(factors):
        m = factors[d]
        parts.append(f"P{d}" + (f"^{m}" if m > 1 else ""))
    return "*".join(parts)


def order_vanishing_at_minus1(poly):
    """Order of vanishing of poly at q=-1 (multiplicity of factor (q+1)=Phi_2)."""
    p = sp.Poly(sp.expand(poly), q)
    if p.is_zero:
        return None
    k = 0
    phi2 = sp.Poly(q + 1, q)
    while True:
        quo, rem = sp.div(p, phi2)
        if rem.is_zero:
            k += 1
            p = quo
        else:
            break
    return k
