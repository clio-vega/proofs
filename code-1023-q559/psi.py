"""Core instrument for Q559: Psi_{(1^r),theta}(M) in exact rational arithmetic.

Two INDEPENDENT implementations:
  psi_headsum  -- directly from eq:Psidef (sum over ordered distinct r-tuples)
  psi_alt      -- via the alternant expansion thm:alt (sum over r-subsets)
Agreement of the two re-verifies thm:alt, on which the whole argument rests.
"""
from fractions import Fraction as F
from itertools import permutations, combinations


def psi_headsum(w, theta):
    """w: list of Fractions (the values of M, in index order). theta: list of ints."""
    r = len(theta)
    n = len(w)
    total = F(0)
    for h in permutations(range(n), r):
        X = [w[i] for i in h]
        term = F(1)
        for j in range(r):
            term *= (1 - 1/X[j]**2) * X[j]**theta[j]
        for j in range(r):
            for k in range(j+1, r):
                term *= (1 - 1/(X[j]*X[k])) / (1 - X[j]/X[k])
        rest = [t for t in range(n) if t not in h]
        for t in rest:
            for j in range(r):
                term *= (1 - 1/(w[t]*X[j])) / (1 - w[t]/X[j])
        total += term
    return total


def B(w, T):
    """B(T) from thm:alt.  T: sorted tuple of indices into w."""
    n = len(w)
    x = [w[i] for i in T]
    out = F(1)
    for xi in x:
        out *= (1 - 1/xi**2)
    for a in range(len(x)):
        for b in range(a+1, len(x)):
            out *= (1 - 1/(x[a]*x[b]))
    for t in range(n):
        if t in T:
            continue
        for xi in x:
            out *= (1 - 1/(w[t]*xi)) / (1 - w[t]/xi)
    return out


def alternant_ratio(x, m):
    """det(x_i^{m_j}) / prod_{i<k}(x_k - x_i)."""
    r = len(m)
    assert len(x) == r
    num = F(0)
    for pi in permutations(range(r)):
        sgn = 1
        # sign of permutation
        seen = [False]*r
        for s in range(r):
            if seen[s]:
                continue
            L = 0
            j = s
            while not seen[j]:
                seen[j] = True
                j = pi[j]
                L += 1
            if L % 2 == 0:
                sgn = -sgn
        p = F(sgn)
        for i in range(r):
            p *= x[i]**m[pi[i]]
        num += p
    den = F(1)
    for i in range(r):
        for k in range(i+1, r):
            den *= (x[k] - x[i])
    return num/den


def psi_alt(w, theta):
    r = len(theta)
    n = len(w)
    m = [theta[j] + j for j in range(r)]   # m_j = theta_j + j - 1, 0-indexed j
    total = F(0)
    for T in combinations(range(n), r):
        x = [w[i] for i in T]
        total += B(w, T) * alternant_ratio(x, m)
    return total


def mvec(theta):
    return [theta[j] + j for j in range(len(theta))]
