"""
Verify: the Schur-basis functionals L_1, L_2, L_3 from the empirical probe at (6,3)
lie in the annihilator of W_6 (i.e., they annihilate all P_ell). This shows the
probe's empirical finding is a manifestation of the general theorem.
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, '/home/clio/projects/probes/2026-08-09-bhattacharya-rhoades-wreath-superspace')
from symmfunc import partitions_of, z_lambda, p_to_s, s_to_p
from sympy import S, Matrix, Rational, Integer, expand
from collections import defaultdict


def epsilon(mu):
    return (-1) ** (sum(mu) - len(mu))


def P_ell_schur(n, ell):
    """P_ell in Schur basis: sum_{l(mu)=ell} eps(mu)/z_mu p_mu, converted."""
    p = {mu: Rational(epsilon(mu), z_lambda(mu))
         for mu in partitions_of(n) if len(mu) == ell}
    return p_to_s(p, n)


def evaluate_L_on_schur(L_dict, F_schur):
    """L is a functional {lam: coeff} extracting Schur coefficients. F is {lam: coeff} in Schur basis."""
    s = S.Zero
    for lam, coef in L_dict.items():
        s += coef * F_schur.get(lam, S.Zero)
    return expand(s)


def main():
    n = 6

    # Probe's functionals in Schur basis:
    L1 = {(4, 2): Rational(8, 3), (3, 3): -8, (3, 2, 1): 1}
    L2 = {(5, 1): -1, (4, 2): 5, (3, 3): -10, (3, 1, 1, 1): 1}
    L3 = {(5, 1): 1, (4, 2): Rational(-20, 3), (3, 3): 15, (2, 2, 2): -5,
          (2, 1, 1, 1, 1): 1}

    for name, L in [('L1', L1), ('L2', L2), ('L3', L3)]:
        print(f"\n--- {name} ---")
        for ell in range(1, n + 1):
            P = P_ell_schur(n, ell)
            v = evaluate_L_on_schur(L, P)
            print(f"  {name}(P_{ell}) = {v}")

    # Also verify L1..L3 evaluated on p_3^2:
    print()
    print("--- L_i on target p_3^2 ---")
    p3sq_s = p_to_s({(3, 3): S.One}, n)
    for name, L in [('L1', L1), ('L2', L2), ('L3', L3)]:
        v = evaluate_L_on_schur(L, p3sq_s)
        print(f"  {name}(p_3^2) = {v}")


if __name__ == "__main__":
    main()
