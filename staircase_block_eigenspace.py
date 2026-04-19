"""
Staircase BLOCK eigenspace test.
Computes the ACTUAL staircase product Φ_j = block_j · block_{j-1} · ... · block_1
and checks whether s_j fixes Φ_j(V^H) at the end of each block j.

block_k = (1+s_k)(1+s_{k-1})...(1+s_1)  [rightmost acts first]

This is different from the sequential product (1+s_j)...(1+s_1) because
blocks repeat earlier generators.
"""

from fractions import Fraction
from itertools import permutations
import time
import sys
import math

def swap_pos(sigma, i):
    """Swap positions i and i+1 (0-indexed)."""
    lst = list(sigma)
    lst[i], lst[i+1] = lst[i+1], lst[i]
    return tuple(lst)

def apply_si_to_vec(vec, i):
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
    rows = []
    pivot_at = {}
    for v in vecs:
        if not v:
            continue
        current = dict(v)
        changed = True
        while changed:
            changed = False
            for col in ordering:
                if col not in current or current[col] == 0:
                    continue
                if col in pivot_at:
                    r_idx = pivot_at[col]
                    factor = current[col] / rows[r_idx][col]
                    for k, c in rows[r_idx].items():
                        current[k] = current.get(k, Fraction(0)) - factor * c
                    current = {k: c for k, c in current.items() if c != 0}
                    changed = True
                    break
        if not current:
            continue
        pivot_col = None
        for col in ordering:
            if col in current and current[col] != 0:
                pivot_col = col
                break
        if pivot_col is None:
            continue
        piv_val = current[pivot_col]
        current = {k: c / piv_val for k, c in current.items()}
        for r_idx in range(len(rows)):
            if pivot_col in rows[r_idx] and rows[r_idx][pivot_col] != 0:
                factor = rows[r_idx][pivot_col]
                for k, c in current.items():
                    rows[r_idx][k] = rows[r_idx].get(k, Fraction(0)) - factor * c
                rows[r_idx] = {k: c for k, c in rows[r_idx].items() if c != 0}
        rows.append(current)
        pivot_at[pivot_col] = len(rows) - 1
    return [r for r in rows if r]

def apply_op(basis, i):
    """Apply (I + s_i) to each vector in basis."""
    result = []
    for v in basis:
        sv = apply_si_to_vec(v, i)
        new_v = vec_add(v, sv)
        if not is_zero(new_v):
            result.append(new_v)
    return result

def check_si_fixes(basis, i):
    for v in basis:
        sv = apply_si_to_vec(v, i)
        diff = vec_sub(sv, v)
        if not is_zero(diff):
            return False
    return True

def compute_VH(n, perms_sorted, ordering):
    H_gens = list(range(0, n-1, 2))  # s_1, s_3, s_5, ... (0-indexed: 0, 2, 4, ...)
    current = [{sigma: Fraction(1)} for sigma in perms_sorted]
    for hi in H_gens:
        raw = apply_op(current, hi)
        current = reduce_basis(raw, ordering)
    return current

def run_test(n):
    print(f"\n{'='*70}")
    print(f"n = {n},  |S_n| = {math.factorial(n)}")

    t0 = time.time()
    perms_sorted = sorted(permutations(range(n)))
    ordering = perms_sorted

    VH = compute_VH(n, perms_sorted, ordering)
    dim_VH = len(VH)
    expected = math.factorial(n) // (2 ** (n // 2))
    print(f"dim(V^H) = {dim_VH}  (expected {expected}, match={dim_VH==expected})")

    # Apply staircase BLOCKS: block_j = (1+s_j)(1+s_{j-1})...(1+s_1)
    # Block j processes generators s_j, s_{j-1}, ..., s_1 (1-indexed)
    # = 0-indexed: j-1, j-2, ..., 0

    current_basis = [dict(v) for v in VH]

    max_check = min(n, 8)  # check s_1 through s_{max_check-1}

    header = f"{'block':>6} {'dim':>5}"
    for k in range(1, max_check):
        header += f" {'s_'+str(k):>7}"
    print(header)
    print("-" * len(header))

    for block_j in range(1, n):  # block 1, 2, ..., n-1 (1-indexed)
        # Apply block_j: generators s_{block_j}, s_{block_j-1}, ..., s_1 (1-indexed)
        # = 0-indexed: block_j-1, block_j-2, ..., 0
        for gen in range(block_j - 1, -1, -1):  # block_j-1 down to 0
            raw = apply_op(current_basis, gen)
            current_basis = reduce_basis(raw, ordering)

        dim_image = len(current_basis)
        line = f"  b={block_j:2}  {dim_image:5}"

        for k in range(1, max_check):
            fixed = check_si_fixes(current_basis, k-1)
            line += f" {'YES' if fixed else 'no':>7}"
        print(line)
        sys.stdout.flush()

    t1 = time.time()
    print(f"\n[n={n} in {t1-t0:.2f}s]")

if __name__ == "__main__":
    print("Staircase BLOCK eigenspace test")
    print("Φ_j = block_j · ... · block_1,  block_k = (1+s_k)...(1+s_1)")
    print("Checking which s_k fix Φ_j(V^H) at end of each block j")

    for n in [3, 4, 5, 6]:
        run_test(n)
        sys.stdout.flush()

    print("\n--- n=7 ---")
    run_test(7)
