"""
Test whether the rank of the BR slice-Frobenius span saturates below dim Lambda_n
as r grows. If yes -> null space persists -> functionals exist for the (i)-theorem.
If it grows to full dim -> obstruction cannot be "linear" in the span sense, needs
different formulation (cone-based, or restricted r).
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, '/home/clio/projects/probes/2026-08-09-bhattacharya-rhoades-wreath-superspace')
from probe import frobenius_op_tilde_slice, sign_twist, target_q3_6
from symmfunc import partitions_of, p_to_s, print_s
from sympy import S, Matrix, Rational, Integer, expand


def slices_frob(n, r):
    """Return list of Schur-basis dicts, one per k=0..n."""
    out = []
    for k in range(n + 1):
        raw = frobenius_op_tilde_slice(n, r, k)
        out.append(sign_twist(raw))
    return out


def rank_and_nullspace(slices, n):
    parts = partitions_of(n)
    d = len(parts)
    ncols = len(slices)
    A = Matrix(d, ncols, lambda i, j: slices[j].get(parts[i], S.Zero))
    r = A.rank()
    N = A.T.nullspace()  # nullspace of A^T = left nullspace of A
    return r, d, A, N


def check_target_hits_nullspace(target_vec, nullspace):
    """target_vec: list of Schur coefficients. nullspace: list of Matrix columns."""
    t = Matrix(target_vec)
    for i, v in enumerate(nullspace):
        val = (v.T * t)[0, 0]
        val = expand(val)
        yield i, val, v


def main():
    n = 6
    parts = partitions_of(n)
    print(f"n = {n}, dim Lambda_n = {len(parts)}")
    print(f"Partitions of {n}: {parts}")
    print()
    # Target = p_3^2
    target = p_to_s({(3, 3): S.One}, n)
    target_vec = [target.get(lam, S.Zero) for lam in parts]

    print("=" * 80)
    print("Test 1: For each fixed r >= 2, rank of individual span")
    print("=" * 80)
    for r in [2, 3, 4, 5, 6, 7]:
        sl = slices_frob(n, r)
        rank, dim, A, N = rank_and_nullspace(sl, n)
        print(f"\n  r = {r}: rank = {rank} / dim = {dim} ({len(N)} left-null vectors)")
        for i, val, v in check_target_hits_nullspace(target_vec, N):
            marker = "  ***target IS detected***" if val != 0 else "  target passes"
            v_str = "[" + ", ".join(str(x) for x in v) + "]"
            print(f"    null #{i}: pairing with target = {val}{marker}")

    print()
    print("=" * 80)
    print("Test 2: COMBINED span across r = 2, 3, ..., R and all k")
    print("=" * 80)
    for R in [3, 4, 5, 6, 7]:
        all_slices = []
        for r in range(2, R + 1):
            all_slices.extend(slices_frob(n, r))
        rank, dim, A, N = rank_and_nullspace(all_slices, n)
        print(f"\n  r_max = {R} (total {len(all_slices)} slices): rank = {rank} / dim = {dim} ({len(N)} null)")
        for i, val, v in check_target_hits_nullspace(target_vec, N):
            marker = "  ***target IS detected***" if val != 0 else "  target passes"
            v_str = "[" + ", ".join(str(x) for x in v) + "]"
            print(f"    null #{i}: pairing = {val}{marker}")
            # Print v as a functional on Schur basis:
            terms = []
            for j, coef in enumerate(v):
                if coef != 0:
                    terms.append(f"{coef} * [s_{parts[j]}]")
            print(f"       functional = " + " + ".join(terms) if terms else "       zero")


if __name__ == "__main__":
    main()
