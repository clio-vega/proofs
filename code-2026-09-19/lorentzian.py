"""Lorentzian-polynomial checker, following Huh-Matherne-Meszaros-St.Dizier
`1906.09633` Definition 2.1 (= Branden-Huh `1902.03719`).

A degree-d homogeneous h = sum_alpha c_alpha x^alpha in n variables is LORENTZIAN iff
  (L1) all c_alpha >= 0,
  (L2) supp(h) is M-convex, and
  (L3) for every beta in N^n with |beta| = d-2, the quadratic form d^beta h has at
       most one positive eigenvalue.

Hessian of d^beta h:   M_ij = c_{beta+e_i+e_j} * (beta+e_i+e_j)!

UNTUNED POSITIVE CONTROL (HMMS Theorem 3): N(s_lambda(x_1..x_m)) is Lorentzian for
every partition lambda.  The checker must pass this before it is aimed at anything new.
"""
from itertools import combinations_with_replacement
from math import factorial
from fractions import Fraction
import numpy as np


# ---------- combinatorics ----------------------------------------------------
def kostka(lam, mu):
    """K_{lam,mu} = # SSYT of shape lam, content mu."""
    lam = tuple(x for x in lam if x > 0)
    if sum(lam) != sum(mu):
        return 0
    m = len(mu)
    cells = [(i, j) for i in range(len(lam)) for j in range(lam[i])]
    filling = {}
    cnt = 0

    def rec(idx, content):
        nonlocal cnt
        if idx == len(cells):
            cnt += 1
            return
        i, j = cells[idx]
        lo = 1
        if j > 0:
            lo = max(lo, filling[(i, j - 1)])
        if i > 0:
            lo = max(lo, filling[(i - 1, j)] + 1)
        for v in range(lo, m + 1):
            if content[v - 1] == mu[v - 1]:
                continue
            filling[(i, j)] = v
            content[v - 1] += 1
            rec(idx + 1, content)
            content[v - 1] -= 1
            del filling[(i, j)]

    rec(0, [0] * m)
    return cnt


def compositions(total, n):
    """All alpha in N^n with |alpha| = total."""
    if n == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, n - 1):
            yield (first,) + rest


def multi_factorial(alpha):
    r = 1
    for a in alpha:
        r *= factorial(a)
    return r


# ---------- M-convexity ------------------------------------------------------
def is_M_convex(support):
    """Symmetric exchange: for alpha,beta in J and i with alpha_i>beta_i, there is j
    with alpha_j<beta_j, alpha-e_i+e_j in J and beta-e_j+e_i in J."""
    J = set(support)
    if not J:
        return True
    n = len(next(iter(J)))
    # all elements must have equal degree for a homogeneous polynomial
    degs = {sum(a) for a in J}
    if len(degs) != 1:
        return False
    for alpha in J:
        for beta in J:
            for i in range(n):
                if alpha[i] <= beta[i]:
                    continue
                ok = False
                for j in range(n):
                    if alpha[j] >= beta[j]:
                        continue
                    a2 = list(alpha); a2[i] -= 1; a2[j] += 1
                    b2 = list(beta);  b2[j] -= 1; b2[i] += 1
                    if tuple(a2) in J and tuple(b2) in J:
                        ok = True
                        break
                if not ok:
                    return False
    return True


# ---------- the checker ------------------------------------------------------
def n_positive_eigs(M, tol=1e-9):
    w = np.linalg.eigvalsh(np.array(M, dtype=float))
    scale = max(1.0, float(np.max(np.abs(w))) if w.size else 1.0)
    return int(np.sum(w > tol * scale))


def is_lorentzian(coeffs, n, verbose=False):
    """coeffs: dict alpha(tuple len n) -> Fraction/float coefficient of x^alpha.
    Returns (bool, reason)."""
    coeffs = {a: c for a, c in coeffs.items() if c != 0}
    if not coeffs:
        return True, "zero polynomial"
    degs = {sum(a) for a in coeffs}
    if len(degs) != 1:
        return False, "not homogeneous"
    d = degs.pop()
    if any(c < 0 for c in coeffs.values()):
        return False, "(L1) negative coefficient"
    if not is_M_convex(set(coeffs)):
        return False, "(L2) support is not M-convex"
    if d < 2:
        return True, "degree < 2: (L3) vacuous"
    for beta in compositions(d - 2, n):
        M = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                g = list(beta); g[i] += 1; g[j] += 1
                g = tuple(g)
                c = coeffs.get(g, 0)
                M[i][j] = float(c) * multi_factorial(g)
        k = n_positive_eigs(M)
        if k > 1:
            if verbose:
                print("   beta =", beta, "has", k, "positive eigenvalues")
            return False, "(L3) d^beta h has %d positive eigenvalues at beta=%s" % (k, beta)
    return True, "Lorentzian"


def normalized_schur(lam, m):
    """N(s_lambda(x_1..x_m)) as a coefficient dict: alpha -> K_{lam,alpha}/alpha!."""
    d = sum(lam)
    out = {}
    for alpha in compositions(d, m):
        K = kostka(lam, alpha)
        if K:
            out[alpha] = Fraction(K, multi_factorial(alpha))
    return out


if __name__ == "__main__":
    print("=== UNTUNED POSITIVE CONTROL: HMMS Thm 3, N(s_lambda) is Lorentzian ===")
    control = [((1,),2),((2,),2),((1,1),2),((2,1),2),((2,1),3),((3,1),3),((2,2),3),
               ((2,2,1),3),((3,2,1),3),((2,1),4),((3,2),3),((4,2,1),3),((2,2,2),3),
               ((3,3),3),((3,1,1),3),((2,2,1,1),4)]
    allok = True
    for lam, m in control:
        h = normalized_schur(lam, m)
        ok, why = is_lorentzian(h, m)
        if not ok:
            allok = False
        print("  %-12s m=%d : %-5s  %s" % (str(lam), m, "PASS" if ok else "FAIL", why))
    print("\nCONTROL PASSED (all normalized Schur Lorentzian):", allok)

    print("\n=== NEGATIVE CONTROL: the detector must be able to say NO ===")
    # x^2+y^2 normalized: not Lorentzian (support {(2,0),(0,2)} is not M-convex)
    bad = {(2,0): Fraction(1,2), (0,2): Fraction(1,2)}
    print("  N(x^2+y^2)          :", is_lorentzian(bad, 2))
    # a bivariate with a log-concavity violation: coeffs (1,0,1) has internal zero
    bad2 = {(3,0):Fraction(1),(2,1):Fraction(1),(1,2):Fraction(0),(0,3):Fraction(1)}
    print("  internal-zero cubic :", is_lorentzian(bad2, 2))
    # HMMS bivariate criterion: a_k^2/C(d,k)^2 >= ...; break it
    bad3 = {(2,0):Fraction(1),(1,1):Fraction(1,100),(0,2):Fraction(1)}
    print("  broken bivariate    :", is_lorentzian(bad3, 2))
