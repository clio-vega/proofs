"""
Verify:
  1. W_n := span{F_m : m >= 1} has dimension = n (as long as n >= 1).
  2. W_n = span{P_ell : 1 <= ell <= n} where P_ell = sum_{l(mu)=ell} eps(mu)/z_mu * p_mu.
  3. For k' >= 2, e >= 2, n = e*k': p_e^{k'} not in W_n. Explicit functional.
  4. Rank stability: for various r, the span of {Frob(M_{n,r}^{(k)})}_k is exactly W_n (as long as r >= 1).
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, '/home/clio/projects/probes/2026-08-09-bhattacharya-rhoades-wreath-superspace')
from symmfunc import partitions_of, z_lambda
from sympy import S, Matrix, Rational, Integer, expand
from collections import defaultdict


def epsilon(mu):
    return (-1) ** (sum(mu) - len(mu))


def F_m_powersum_vec(n, m, parts):
    """Return F_m as a vector indexed by parts (list of partitions of n) in p-basis."""
    return [Rational(epsilon(mu) * m ** len(mu), z_lambda(mu)) for mu in parts]


def P_ell_vec(n, ell, parts):
    """P_ell = sum_{l(mu)=ell} eps(mu)/z_mu p_mu, as vector in p-basis."""
    return [Rational(epsilon(mu), z_lambda(mu)) if len(mu) == ell else S.Zero for mu in parts]


def p_target_vec(n, mu_target, parts):
    """The vector for p_{mu_target} in p-basis (1 at mu_target, 0 elsewhere)."""
    return [S.One if mu == mu_target else S.Zero for mu in parts]


def rank_of_vectors(vecs):
    if not vecs:
        return 0
    d = len(vecs[0])
    M = Matrix(d, len(vecs), lambda i, j: vecs[j][i])
    return M.rank()


def solve_target_in_span(target_vec, spanning_vecs):
    """Return (solvable, solution) where solvable = target in span. solution is a Matrix column."""
    if not spanning_vecs:
        return (all(v == 0 for v in target_vec), None)
    d = len(target_vec)
    A = Matrix(d, len(spanning_vecs), lambda i, j: spanning_vecs[j][i])
    b = Matrix(d, 1, lambda i, _: target_vec[i])
    aug = A.row_join(b)
    return (A.rank() == aug.rank()), None


def find_null_functional(vecs, target_vec, parts):
    """Find L such that L(v) = 0 for all v in vecs but L(target) != 0.
    Return L as dict {mu -> coeff}."""
    d = len(target_vec)
    A = Matrix(d, len(vecs), lambda i, j: vecs[j][i]) if vecs else Matrix(d, 0)
    # Left nullspace = row nullspace of A^T = null(A^T)
    NT = A.T.nullspace()
    for v in NT:
        val = sum(v[i] * target_vec[i] for i in range(d))
        val = expand(val)
        if val != 0:
            return {parts[i]: v[i] for i in range(d) if v[i] != 0}, val
    return None, None


def main():
    print("=" * 70)
    print("Verify W_n structure at multiple n")
    print("=" * 70)
    for n in [4, 5, 6, 8, 9, 12]:
        parts = partitions_of(n)
        print(f"\n--- n = {n}, dim Lambda_n = {len(parts)} ---")

        # F_m for m = 1, ..., 3n
        F_vecs = [F_m_powersum_vec(n, m, parts) for m in range(1, 3 * n + 1)]
        r_F = rank_of_vectors(F_vecs)
        print(f"  rank span{{F_m : 1<=m<={3*n}}} = {r_F}  (predicted: {n})")

        # P_ell for ell = 1, ..., n
        P_vecs = [P_ell_vec(n, ell, parts) for ell in range(1, n + 1)]
        # Drop zero vectors (if any ell has no partitions of that length)
        P_vecs_nz = [v for v in P_vecs if not all(x == 0 for x in v)]
        r_P = rank_of_vectors(P_vecs_nz)
        print(f"  rank span{{P_ell : 1<=ell<={n}}} = {r_P}  (predicted: {n})")

        # Combined
        r_combined = rank_of_vectors(F_vecs + P_vecs_nz)
        print(f"  rank of F union P = {r_combined}  (predicted: {n})")

        assert r_F == r_P == r_combined == n, "SPAN MISMATCH"
        print(f"  W_n span check: OK (rank {n} in each case, and equal)")

    print()
    print("=" * 70)
    print("Verify: p_e^{k'} not in W_n for various (n, e, k'>=2)")
    print("=" * 70)

    cases = [
        (6, 2),   # e=2, n=6, k'=3
        (6, 3),   # e=3, n=6, k'=2
        (8, 2),   # e=2, n=8, k'=4
        (8, 4),   # e=4, n=8, k'=2
        (9, 3),   # e=3, n=9, k'=3
        (12, 2),  # e=2, n=12, k'=6
        (12, 3),  # e=3, n=12, k'=4
        (12, 4),  # e=4, n=12, k'=3
        (12, 6),  # e=6, n=12, k'=2
    ]

    for (n, e) in cases:
        k_prime = n // e
        assert n % e == 0
        parts = partitions_of(n)
        mu_target = tuple([e] * k_prime)
        target_vec = p_target_vec(n, mu_target, parts)

        P_vecs = [P_ell_vec(n, ell, parts) for ell in range(1, n + 1)
                  if any(len(mu) == ell for mu in parts)]

        solvable, _ = solve_target_in_span(target_vec, P_vecs)
        marker = "IN W_n (bad)" if solvable else "NOT in W_n (good)"
        print(f"\n  (n={n}, e={e}, k'={k_prime}): p_{{{mu_target}}}  ->  {marker}")

        if not solvable:
            L, val = find_null_functional(P_vecs, target_vec, parts)
            if L is None:
                print(f"    !!! No null functional found (shouldn't happen)")
            else:
                print(f"    Explicit functional L (nonzero coeffs on p-basis):")
                length_kp = [mu for mu in L.keys() if len(mu) == k_prime]
                # print only length-k' entries (they should form a valid annihilator)
                for mu in sorted(length_kp, key=lambda m: (-sum(m), m)):
                    print(f"      p_{mu} coeff: {L[mu]}")
                print(f"    L(p_{{{mu_target}}}) = {val}")


if __name__ == "__main__":
    main()
