"""
Generalized d-core / d-quotient via abacus (d runners), plus reuse of
SYT/G_lambda/character machinery from the existing 2-quotient library.
SageMath NOT available -> everything self-contained with sympy.
"""
import sys
sys.path.insert(0, '/home/clio/projects/scratch/2026-06-04-sT-2quotient')
from functools import lru_cache
import sympy as sp
from sympy import symbols, binomial, factorial, Rational

from syt_lib import (partitions, syt_list, descent_set, maj, s_stat,
                     G_lambda, fake_degree, num_syt, n_of_lambda,
                     two_core_and_quotient, beta_set)

q = symbols('q')


# --- Murnaghan-Nakayama character chi^lam(mu) ----------------------------
@lru_cache(maxsize=None)
def mn_char(lam, mu):
    lam = tuple(x for x in lam if x > 0)
    if sum(lam) == 0 and len(mu) == 0:
        return 1
    if len(mu) == 0:
        return 1 if sum(lam) == 0 else 0
    r = mu[0]; rest = mu[1:]; total = 0; L = len(lam)
    if L == 0:
        return 0
    beta = [lam[i] + (L - 1 - i) for i in range(L)]; bset = set(beta)
    for b in beta:
        if b - r >= 0 and (b - r) not in bset:
            nb = sorted([x for x in beta if x != b] + [b - r], reverse=True)
            nl = tuple(nb[i] - (L - 1 - i) for i in range(L))
            nl = tuple(x for x in nl if x > 0)
            h = sum(1 for x in beta if b - r < x < b)
            total += ((-1) ** h) * mn_char(nl, rest)
    return total


# --- d-core and d-quotient via abacus (d runners) ------------------------
def d_core_and_quotient(lam, d):
    """
    Returns (core_partition, (q0,...,q_{d-1})) for any d>=2.
    Convention: beta_i = lam_i + (r-1-i), r beads with r divisible by d,
    runner = beta mod d, row = beta // d.  Standard James/Kerber.
    """
    lam = list(lam)
    r = len(lam)
    if r % d != 0:
        r += d - (r % d)
    r = max(r, d)
    beta = beta_set(lam, r)  # length r, strictly decreasing

    runners = {k: [] for k in range(d)}
    for b in beta:
        runners[b % d].append(b // d)
    for k in runners:
        runners[k].sort(reverse=True)

    quotient = []
    for k in range(d):
        vals = runners[k]
        m = len(vals)
        part = [vals[i] - (m - 1 - i) for i in range(m)]
        part = tuple(x for x in part if x > 0)
        quotient.append(part)

    # core: push beads up on each runner
    core_beta = []
    for k in range(d):
        m = len(runners[k])
        for row in range(m):
            core_beta.append(k + d * row)
    core_beta.sort(reverse=True)
    rc = len(core_beta)
    core = [core_beta[i] - (rc - 1 - i) for i in range(rc)]
    core = tuple(x for x in core if x > 0)
    return core, tuple(quotient)


def multinom(parts):
    tot = sum(parts)
    res = factorial(tot)
    for p in parts:
        res /= factorial(p)
    return int(res)


def G_minus1(lam):
    """G_lambda(-1) as an integer."""
    return int(sp.Poly(G_lambda(lam), q).eval(-1))


def ord_at_root(poly, zeta_minimal_poly):
    """Order of vanishing of poly at a root of zeta_minimal_poly (irreducible)."""
    p = sp.Poly(sp.expand(poly), q)
    if p.is_zero:
        return None
    phi = sp.Poly(zeta_minimal_poly, q)
    k = 0
    while True:
        quo, rem = sp.div(p, phi)
        if rem.is_zero:
            k += 1
            p = quo
        else:
            break
    return k


if __name__ == "__main__":
    # sanity: d=2 matches existing routine
    bad = 0
    for n in range(1, 11):
        for lam in partitions(n):
            c2, q2 = d_core_and_quotient(lam, 2)
            c2b, q2b = two_core_and_quotient(lam)
            if c2 != c2b or set(q2) != set(q2b):
                # quotient ordering may differ; compare as multiset of partitions
                if sorted(q2) != sorted(q2b) or c2 != c2b:
                    bad += 1
                    print("MISMATCH d=2", lam, (c2, q2), (c2b, q2b))
    print("d=2 self-check mismatches:", bad)
