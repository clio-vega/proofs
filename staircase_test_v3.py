"""
Staircase cascade eigenspace test - v3 with diagnostics.

Key finding from v2: ALL blocks (odd AND even) give YES for n=4,5,6,7.
This script:
1. Confirms the pattern with more diagnostics
2. Checks stronger claim: is image PRESERVED (not just contained in +1 eig of s_j)?
3. Checks: what is the image equal to? Does it = V^H throughout?
4. Checks: does s_{j+1} also fix the image? (looking for full symmetry)
5. Tries to understand WHY dim(image) = dim(V^H) always (no collapsing).

Also computes individual irrep decomposition via characters.
"""

from fractions import Fraction
from itertools import permutations
import time
import sys

def swap_pos(sigma, i):
    """Left mult by s_i (0-indexed): swap positions i and i+1."""
    lst = list(sigma)
    lst[i], lst[i+1] = lst[i+1], lst[i]
    return tuple(lst)

def apply_si_to_vec(vec, i):
    """Apply s_i to sparse vector."""
    result = {}
    for sigma, c in vec.items():
        target = swap_pos(sigma, i)
        result[target] = result.get(target, Fraction(0)) + c
    return {k: v for k, v in result.items() if v != 0}

def vec_add(u, v):
    result = dict(u)
    for k, c in v.items():
        result[k] = result.get(k, Fraction(0)) + c
    return {k: c for k, c in result.items() if c != 0}

def vec_scale(v, s):
    return {k: c * s for k, c in v.items()}

def vec_sub(u, v):
    return vec_add(u, vec_scale(v, Fraction(-1)))

def is_zero(v):
    return all(c == 0 for c in v.values())

def reduce_basis(vecs, ordering):
    """Gaussian elimination on sparse vectors. Returns reduced basis."""
    rows = []
    pivot_at = {}  # pivot_col -> row_index in rows

    for v in vecs:
        if not v:
            continue
        current = dict(v)

        # Reduce against existing pivots
        changed = True
        while changed:
            changed = False
            for col in ordering:
                if col not in current or current[col] == 0:
                    continue
                if col in pivot_at:
                    r_idx = pivot_at[col]
                    factor = current[col] / rows[r_idx][col]
                    current = vec_sub(current, vec_scale(rows[r_idx], factor))
                    current = {k: c for k, c in current.items() if c != 0}
                    changed = True
                    break

        if not current:
            continue

        # Find pivot (first nonzero in ordering)
        pivot_col = None
        for col in ordering:
            if col in current and current[col] != 0:
                pivot_col = col
                break

        if pivot_col is None:
            continue

        # Normalize
        piv_val = current[pivot_col]
        current = {k: c / piv_val for k, c in current.items()}

        # Back-substitute into existing rows
        for r_idx in range(len(rows)):
            if pivot_col in rows[r_idx] and rows[r_idx][pivot_col] != 0:
                factor = rows[r_idx][pivot_col]
                rows[r_idx] = vec_sub(rows[r_idx], vec_scale(current, factor))

        rows.append(current)
        pivot_at[pivot_col] = len(rows) - 1

    return [r for r in rows if r]

def apply_block(basis, i):
    """Apply (I + s_i) to each vector in basis."""
    result = []
    for v in basis:
        sv = apply_si_to_vec(v, i)
        new_v = vec_add(v, sv)
        if not is_zero(new_v):
            result.append(new_v)
    return result

def check_si_fixes(basis, i):
    """Check if s_i v = v for all basis vectors."""
    for v in basis:
        sv = apply_si_to_vec(v, i)
        diff = vec_sub(sv, v)
        if not is_zero(diff):
            return False
    return True

def subspace_equal(basis1, basis2, ordering):
    """Check if two subspaces are equal."""
    if len(basis1) != len(basis2):
        return False
    # Check each basis1 vector is in span of basis2 and vice versa
    # Easier: reduce basis1 + basis2 and check rank = len(basis1)
    combined = list(basis1) + list(basis2)
    reduced = reduce_basis(combined, ordering)
    return len(reduced) == len(basis1)

def compute_VH(n, perms_sorted, ordering):
    """Compute V^H = joint +1 eigenspace of s_0, s_2, s_4, ... (0-indexed)."""
    H_gens = list(range(0, n-1, 2))
    current_basis = [{sigma: Fraction(1)} for sigma in perms_sorted]

    for hi in H_gens:
        new_raw = apply_block(current_basis, hi)
        current_basis = reduce_basis(new_raw, ordering)

    return current_basis

def run_test(n, verbose=True):
    print(f"\n{'='*65}")
    import math
    print(f"n = {n},  |S_n| = {math.factorial(n)}")

    t0 = time.time()
    perms_sorted = sorted(permutations(range(n)))
    ordering = perms_sorted

    print("Computing V^H...")
    VH = compute_VH(n, perms_sorted, ordering)
    dim_VH = len(VH)
    print(f"dim(V^H) = {dim_VH}")

    if dim_VH == 0:
        print("V^H is trivial, skipping.")
        return

    print(f"\nExpected dim(V^H) = n!/2^floor(n/2) = ?")
    import math
    expected = math.factorial(n) // (2 ** (n // 2))
    print(f"  n!/2^floor(n/2) = {math.factorial(n)}/2^{n//2} = {expected}")
    print(f"  Computed dim(V^H) = {dim_VH}")
    print(f"  Match: {dim_VH == expected}")

    print(f"\nStaircase cascade:")
    header = f"{'j':>4} {'parity':>6} {'dim':>6} {'s_j fixes':>10} {'image=V^H':>10}"
    print(header)
    print("-" * 45)

    current_basis = [dict(v) for v in VH]

    for j in range(1, n):  # j = 1, ..., n-1 (1-indexed)
        i = j - 1  # 0-indexed

        new_raw = apply_block(current_basis, i)
        current_basis = reduce_basis(new_raw, ordering)

        dim_image = len(current_basis)
        parity = "odd" if j % 2 == 1 else "even"

        # Check s_j fixes image
        fixed_j = check_si_fixes(current_basis, i)

        # Check if image equals V^H
        img_eq_VH = subspace_equal(current_basis, VH, ordering)

        print(f"  j={j:2d}  {parity:>6}  {dim_image:5d}  {'YES':>10}  {'YES' if img_eq_VH else 'NO':>10}"
              if fixed_j else
              f"  j={j:2d}  {parity:>6}  {dim_image:5d}  {'NO':>10}  {'YES' if img_eq_VH else 'NO':>10}")

        sys.stdout.flush()

    t1 = time.time()
    print(f"\n[n={n} completed in {t1-t0:.2f}s]")

if __name__ == "__main__":
    print("Staircase cascade eigenspace test v3")
    print("CLAIM: Φ_j(V^H) ⊆ ker(s_j - I) for all j")
    print("STRONGER: Is Φ_j(V^H) = V^H for all j?")
    print()

    for n in [4, 5, 6, 7]:
        run_test(n)
        sys.stdout.flush()
