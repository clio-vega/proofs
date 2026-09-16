"""
PROVE 2026-08-11: Verify two structural claims empirically.

Claim (A): the coefficient of |e*lambda> in v_{k',e} mod qL is f^lambda
           (number of standard Young tableaux of shape lambda), and the
           support consists ONLY of partitions of the form e*lambda for
           lambda \\vdash k'.

Claim (B) (INDIVIDUAL): for every partition lambda and every i in Z/eZ,
          epsilon_i(e*lambda) <= 1 (individual crystal element bound).

Broader sweep than the 12 target triples in epsilon_i_probe.py.
"""

import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "2026-08-12-q34-jacon-lacabanne"))

from qfock import v_kprime_e, vec_support, partitions
from probe4_canonical import good_remove_i_node, remove_node
from epsilon_i_probe import v_at_q0, tilde_e_i_at_q0, epsilon_i


def f_lambda(lam):
    """Number of standard Young tableaux of shape lam, via hook-length formula."""
    if not lam:
        return 1
    n = sum(lam)
    # hooks
    prod = 1
    L = len(lam)
    for r in range(L):
        for c in range(lam[r]):
            # arm: cells strictly right in row r
            arm = lam[r] - 1 - c
            # leg: cells strictly below in column c
            leg = 0
            for rr in range(r + 1, L):
                if lam[rr] > c:
                    leg += 1
                else:
                    break
            prod *= (arm + leg + 1)
    from math import factorial
    return factorial(n) // prod


def e_scale(lam, e):
    return tuple(e * x for x in lam)


def check_expansion(k_prime, e):
    """Verify (A) for a given (k', e)."""
    v = v_kprime_e(k_prime, e)
    v0 = v_at_q0(v)
    # Expected: v0[e*lambda] = f^lambda for each lambda vdash k'.
    #           No other partitions in support.
    ok = True
    expected_support = set()
    for lam in partitions(k_prime):
        mu = e_scale(lam, e)
        expected_support.add(mu)
        expected_coeff = f_lambda(lam)
        actual_coeff = v0.get(mu, 0)
        if actual_coeff != expected_coeff:
            print(f"  MISMATCH: e={e}, k'={k_prime}, lambda={lam}, "
                  f"e*lambda={mu}: expected f^lambda={expected_coeff}, "
                  f"got {actual_coeff}")
            ok = False
    # Check no extra support
    for mu, c in v0.items():
        if mu not in expected_support and c != 0:
            print(f"  EXTRA SUPPORT: e={e}, k'={k_prime}, mu={mu} with coeff {c}")
            ok = False
    return ok


def check_individual(lam, e, max_cap=5):
    """Check epsilon_i(e*lambda) <= 1 for all i in Z/eZ."""
    mu = e_scale(lam, e)
    all_ok = True
    max_eps = 0
    for i in range(e):
        v0 = {mu: 1}
        eps, _ = epsilon_i(v0, i, e, cap=max_cap)
        if eps > 1:
            print(f"  MISMATCH: e={e}, lambda={lam} -> e*lambda={mu}, i={i}: "
                  f"epsilon_i = {eps} > 1")
            all_ok = False
        max_eps = max(max_eps, eps)
    return all_ok, max_eps


def main():
    print("=" * 78)
    print("(A) Expansion check: [v_{k',e}] mod qL = sum_{lam vdash k'} f^lam |e*lam>")
    print("=" * 78)
    total_A = 0
    pass_A = 0
    for k_prime in range(1, 6):
        for e in range(2, 6):
            if e * k_prime > 16:
                continue
            total_A += 1
            ok = check_expansion(k_prime, e)
            status = "OK" if ok else "FAIL"
            print(f"  (e, k') = ({e}, {k_prime}), |lambda| = {e*k_prime}: {status}")
            if ok:
                pass_A += 1
    print(f"\n(A) TOTAL: {pass_A}/{total_A} PASS")

    print()
    print("=" * 78)
    print("(B) Individual claim: epsilon_i(e*lambda) <= 1 for all i, all lambda")
    print("=" * 78)
    total_B = 0
    pass_B = 0
    max_eps_seen = 0
    counts = {0: 0, 1: 0}
    for k_prime in range(1, 8):
        for e in range(2, 6):
            if e * k_prime > 20:
                continue
            for lam in partitions(k_prime):
                total_B += 1
                ok, max_eps = check_individual(lam, e)
                if ok:
                    pass_B += 1
                    max_eps_seen = max(max_eps_seen, max_eps)
                    for i in range(e):
                        v0 = {e_scale(lam, e): 1}
                        eps, _ = epsilon_i(v0, i, e)
                        counts[eps] = counts.get(eps, 0) + 1
    print(f"\n(B) TOTAL: {pass_B}/{total_B} PASS (max epsilon_i seen: {max_eps_seen})")
    print(f"    Distribution of epsilon_i values: {dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
