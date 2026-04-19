"""
Staircase cascade eigenspace test - optimized version.

Key insight: work with the PERMUTATION GROUP ALGEBRA using sparse representations.
Instead of N×N matrices, represent elements of the group algebra as dicts {perm: coeff}.

s_i acts on the LEFT: s_i · (sum c_σ e_σ) = sum c_σ e_{s_i σ}
So the action on basis elements e_σ sends e_σ → e_{s_i σ}.

V^H = joint +1 eigenspace of s_1, s_3, ... (1-indexed)

For the staircase: we don't need full N×N matrices.
We work with SUBSPACES represented as spans of group algebra elements.

But for exact arithmetic, let's use a smarter representation:
- Represent each basis vector as a dict {perm_tuple: Fraction}
- Gaussian elimination on these sparse vectors

This avoids building N×N matrices explicitly.
"""

from fractions import Fraction
from itertools import permutations
import time
import sys

def swap_pos(sigma, i):
    """Left multiplication by s_i: swap positions i and i+1 (0-indexed). Returns new tuple."""
    lst = list(sigma)
    lst[i], lst[i+1] = lst[i+1], lst[i]
    return tuple(lst)

def apply_si_to_vec(vec, i):
    """Apply s_i (0-indexed) to a sparse vector (dict perm->Fraction).
    s_i acts by left multiplication: e_sigma -> e_{s_i sigma}.
    """
    result = {}
    for sigma, c in vec.items():
        target = swap_pos(sigma, i)
        result[target] = result.get(target, Fraction(0)) + c
    return {k: v for k, v in result.items() if v != 0}

def vec_add(u, v):
    result = dict(u)
    for k, c in v.items():
        result[k] = result.get(k, Fraction(0)) + c
    return {k: v for k, v in result.items() if v != 0}

def vec_scale(v, s):
    return {k: c * s for k, c in v.items()}

def vec_sub(u, v):
    return vec_add(u, vec_scale(v, Fraction(-1)))

def vec_dot_norm(v):
    """Check if vector is zero."""
    return all(c == 0 for c in v.values())

def apply_block(basis, i):
    """Apply (I + s_i) to each vector in basis, return new list of vectors."""
    new_basis = []
    for v in basis:
        sv = apply_si_to_vec(v, i)
        new_v = vec_add(v, sv)  # (I + s_i) v
        if not vec_dot_norm(new_v):
            new_basis.append(new_v)
    return new_basis

def reduce_basis(vecs, ordering):
    """
    Gaussian elimination on sparse vectors.
    ordering: list of perm tuples giving column order (for pivot selection).
    Returns reduced basis.
    """
    # Use a pivot structure: for each row, track which column is pivot
    rows = []
    pivot_at = {}  # col -> row_index

    for v in vecs:
        if not v:
            continue
        # Find first nonzero entry in ordering
        current = dict(v)

        # Reduce against existing pivots
        for col in ordering:
            if col not in current or current[col] == 0:
                continue
            if col in pivot_at:
                row_idx = pivot_at[col]
                factor = current[col] / rows[row_idx][col]
                current = vec_sub(current, vec_scale(rows[row_idx], factor))
                # clean zeros
                current = {k: cv for k, cv in current.items() if cv != 0}

        if not current:
            continue

        # Find pivot of current (first nonzero in ordering)
        pivot_col = None
        for col in ordering:
            if col in current and current[col] != 0:
                pivot_col = col
                break

        if pivot_col is None:
            continue

        # Normalize so pivot = 1
        piv_val = current[pivot_col]
        current = {k: cv / piv_val for k, cv in current.items()}

        # Eliminate pivot_col from existing rows
        for r_idx, row in enumerate(rows):
            if pivot_col in row and row[pivot_col] != 0:
                factor = row[pivot_col]
                rows[r_idx] = {k: cv - factor * current.get(k, Fraction(0))
                               for k, cv in row.items()}
                rows[r_idx].update({k: -factor * current.get(k, Fraction(0))
                                   for k in current if k not in row})
                rows[r_idx] = {k: cv for k, cv in rows[r_idx].items() if cv != 0}

        rows.append(current)
        pivot_at[pivot_col] = len(rows) - 1

    return [r for r in rows if r]

def compute_VH(n, perms_sorted, ordering):
    """
    Compute V^H = joint +1 eigenspace of s_1, s_3, ... (1-indexed = s_0, s_2, ... 0-indexed).

    +1 eigenspace of s_i: vectors v such that s_i v = v.
    A basis: for each orbit {sigma, s_i sigma} with sigma < s_i sigma,
      take e_sigma + e_{s_i sigma} (if orbit size 2)
      or e_sigma (if sigma is fixed by s_i, i.e., s_i sigma = sigma).

    The +1 eigenspace of s_i has basis:
      {e_sigma + e_{s_i sigma} : sigma < s_i sigma} ∪ {e_sigma : s_i sigma = sigma}

    Intersection of these for all H generators.

    We compute intersection by:
    1. Start with the +1 eigenspace of s_0
    2. Intersect with +1 eigenspace of s_2, etc.
    """

    H_gens = list(range(0, n-1, 2))  # 0-indexed: 0, 2, 4, ...

    # Start with full space (basis = all e_sigma)
    # Represent as sparse vectors
    current_basis = [{sigma: Fraction(1)} for sigma in perms_sorted]

    for hi in H_gens:
        # Project onto +1 eigenspace of s_hi
        # P_{+1} = (I + s_hi) / 2
        # Apply (I + s_hi) to each basis vector, then re-reduce
        new_raw = []
        for v in current_basis:
            sv = apply_si_to_vec(v, hi)
            proj = vec_add(v, sv)  # (I + s_hi) v, unnormalized (scale by 2)
            if not vec_dot_norm(proj):
                new_raw.append(proj)

        current_basis = reduce_basis(new_raw, ordering)
        print(f"  After projecting onto +1 eig of s_{hi+1}: dim = {len(current_basis)}")

    return current_basis

def check_si_fixes(basis, i):
    """Check if s_i fixes every vector in basis."""
    for v in basis:
        sv = apply_si_to_vec(v, i)
        diff = vec_sub(sv, v)
        if not vec_dot_norm(diff):
            return False
    return True

def run_test(n):
    print(f"\n{'='*60}")
    import math
    fact_n = math.factorial(n)
    print(f"n = {n},  |S_n| = {fact_n}")

    t0 = time.time()
    perms_sorted = sorted(permutations(range(n)))
    ordering = perms_sorted  # Use lex order for pivoting

    print("Computing V^H...")
    VH = compute_VH(n, perms_sorted, ordering)
    print(f"dim(V^H) = {len(VH)}")

    if len(VH) == 0:
        print("V^H is trivial, skipping.")
        return

    print(f"\nStaircase cascade:")
    print(f"{'Block j':>8} {'parity':>8} {'dim(image)':>12} {'s_j fixes image?':>18}")
    print("-" * 55)

    current_basis = [dict(v) for v in VH]

    for j in range(1, n):  # j = 1, ..., n-1 (1-indexed)
        i = j - 1  # 0-indexed transposition

        # Apply (I + s_i) to each vector in current_basis
        new_raw = apply_block(current_basis, i)

        # Reduce to basis
        current_basis = reduce_basis(new_raw, ordering)

        dim_image = len(current_basis)

        # Check: does s_j (= Si[i]) fix every vector?
        fixed = check_si_fixes(current_basis, i)

        parity = "odd" if j % 2 == 1 else "even"
        fixed_str = "YES" if fixed else "NO"

        print(f"  j={j:2d}   {parity:>6}   dim={dim_image:5d}   s_{j} fixes image: {fixed_str}")

    t1 = time.time()
    print(f"[n={n} completed in {t1-t0:.2f}s]")

if __name__ == "__main__":
    print("Staircase cascade eigenspace test (sparse exact arithmetic)")
    print("Claim: After Φ_j applied to V^H, image lies in +1 eigenspace of s_j")
    print()

    for n in [4, 5, 6, 7]:
        run_test(n)
        sys.stdout.flush()

    print("\nAttempting n=8...")
    run_test(8)
