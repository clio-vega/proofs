"""
PROVE 2026-08-11 Phase 2: empirically compute epsilon_i(v_{k',e}) for the 12
target triples.

Strategy:
    1. Compute v = v_{k',e} in the standard basis of level-1 Uglov q-Fock
       (using qfock.v_kprime_e). All coefficients are in Z[q] (nonneg powers).
    2. Reduce mod qL: replace each coefficient by its q=0 constant term.
       This gives a Z-linear combination of standard basis vectors, which
       is exactly the class of v in L/qL (level-1: standard-at-q=0 = crystal-at-q=0).
    3. Iteratively apply Kashiwara tilde_e_i via the good-removable-i-node rule:
         tilde_e_i [mu] = [mu \\ good-i-node] if that node exists, else 0.
       reusing probe4_canonical.good_remove_i_node.
    4. Count how many iterations until zero. That's epsilon_i(v_{k',e}).

Reports: PASS if all 12 triples have epsilon_i <= 1; FAIL otherwise.
"""

import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "2026-08-12-q34-jacon-lacabanne"))
sys.path.insert(0, os.path.join(HERE, "..", "2026-08-13-iijima-B1"))

from qfock import (
    v_kprime_e, vec_support, vec_zero,
)
from probe4_canonical import good_remove_i_node, remove_node


# ---------------------------------------------------------------------------
# Kashiwara-lattice reduction and tilde_e_i at q = 0
# ---------------------------------------------------------------------------

def v_at_q0(v):
    """Reduce Fock vector v mod qL: replace each Laurent-poly coeff by its q=0
    value (the coefficient of q^0). Returns dict {partition: int}.

    Prerequisite: for each mu in supp(v), coeff(mu) must lie in Z[[q]] (no q^-1
    poles), else the vector is not in L and the class is undefined.
    """
    out = {}
    for lam, poly in v.items():
        if not poly:
            continue
        min_deg = min(poly.keys())
        if min_deg < 0:
            raise ValueError(f"Coefficient of |{lam}> has q^{min_deg} term — "
                             f"not in Kashiwara lattice L.")
        c0 = poly.get(0, 0)
        if c0 != 0:
            out[lam] = c0
    return out


def tilde_e_i_at_q0(v0, i, e):
    """Apply Kashiwara tilde_e_i to a q=0 vector.

    v0 : dict {partition: int}    (an element of L/qL, expressed in the
                                    crystal-basis = standard-basis-at-q=0)
    i  : residue in {0, 1, ..., e-1}
    e  : level

    Returns dict {partition: int} — coefficient of each surviving partition
    after applying tilde_e_i.
    """
    out = {}
    for mu, c in v0.items():
        if c == 0:
            continue
        good = good_remove_i_node(mu, i, e)
        if good is None:
            continue  # tilde_e_i annihilates this crystal element
        nu = remove_node(mu, good)
        out[nu] = out.get(nu, 0) + c
        if out[nu] == 0:
            del out[nu]
    return out


def epsilon_i(v0, i, e, cap=10):
    """Compute epsilon_i of a q=0 vector: max n >= 0 such that tilde_e_i^n v0 != 0.

    Also returns the chain (v0, tilde_e_i v0, tilde_e_i^2 v0, ...) up to the
    last nonzero one and then a final 0.
    """
    chain = [dict(v0)]
    cur = v0
    n = 0
    while cur and n < cap:
        cur = tilde_e_i_at_q0(cur, i, e)
        chain.append(cur)
        if not cur:
            return n, chain
        n += 1
    return n, chain


def fmt_v0(v0, max_terms=6):
    """Compact string form for a q=0 vector."""
    if not v0:
        return "0"
    parts = []
    items = sorted(v0.items(), key=lambda kv: (-sum(kv[0]), tuple(-x for x in kv[0])))
    for mu, c in items[:max_terms]:
        sign = "+" if c > 0 else "-"
        abs_c = abs(c)
        if abs_c == 1:
            parts.append(f"{sign}|{mu}>")
        else:
            parts.append(f"{sign}{abs_c}|{mu}>")
    s = " ".join(parts)
    if s.startswith("+"):
        s = s[1:]
    if len(items) > max_terms:
        s += f" ... ({len(items) - max_terms} more)"
    return s


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    tests = [
        (2, 2), (2, 3),
        (3, 2), (3, 3),
        (4, 2),
    ]

    print("=" * 78)
    print("PROVE 2026-08-11 Phase 2:  epsilon_i(v_{k',e}) via good-node rule")
    print("=" * 78)

    total_triples = 0
    pass_triples = 0
    fail_details = []

    for (e, kp) in tests:
        print()
        print(f"### (e, k') = ({e}, {kp}) :  v_{{{kp},{e}}} = P_{e}^{kp} |emptyset>")
        print(f"    (partition weight |lambda| = {e * kp})")

        v = v_kprime_e(kp, e)
        supp_size = len(vec_support(v))
        print(f"    v has standard-basis support of size {supp_size}.")

        try:
            v0 = v_at_q0(v)
        except ValueError as ex:
            print(f"    !!! {ex}")
            print(f"    !!! v is NOT in the Kashiwara lattice L. C5(a) FAILS at this triple.")
            for i in range(e):
                fail_details.append((e, kp, i, "not-in-L"))
            continue

        print(f"    [v] mod qL has support of size {len(v0)}.")
        print(f"    [v] = {fmt_v0(v0)}")

        for i in range(e):
            total_triples += 1
            eps, chain = epsilon_i(v0, i, e)
            status = "PASS" if eps <= 1 else "FAIL"
            marker = "  <-- FAIL" if eps > 1 else ""
            print(f"      i = {i}:  epsilon_{i}(v) = {eps}  [{status}]{marker}")
            for k, c in enumerate(chain[1:], start=1):
                print(f"          tilde_e_{i}^{k} [v] = {fmt_v0(c)}")
            if eps <= 1:
                pass_triples += 1
            else:
                fail_details.append((e, kp, i, f"epsilon={eps}"))

    print()
    print("=" * 78)
    print(f"OVERALL:  {pass_triples}/{total_triples} triples PASS (epsilon_i <= 1)")
    print("=" * 78)
    if fail_details:
        print("Failing triples:")
        for (e, kp, i, why) in fail_details:
            print(f"  (e={e}, k'={kp}, i={i}): {why}")
    else:
        print("*** ALL PASS — C5 empirical verification complete for level 1 target set. ***")


if __name__ == "__main__":
    main()
