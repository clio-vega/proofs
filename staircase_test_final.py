"""
Staircase cascade eigenspace test - FINAL VERSION.
Tests: after Φ_j = block_1 · ... · block_j applied to V^H,
does s_j fix every vector in the resulting image?

SETUP (1-indexed):
- s_i swaps positions i and i+1 in one-line notation (left multiplication)
- block_j = (1+s_j)(1+s_{j-1})...(1+s_1)  [rightmost acts first on input]
  Wait: re-read the spec: block j = (1+s_j)...(1+s_1), Φ_j = block_1 · block_2 · ... · block_j
  So block_1 = (1+s_1), block_2 = (1+s_2)(1+s_1), block_j = (1+s_j)...(1+s_1)
  Φ_j = block_1 · block_2 · ... · block_j

  Applied to a vector v:
  Φ_j(v) = block_j(block_{j-1}(...(block_1(v))...))

  block_1(v) = (1+s_1)v
  block_2(block_1(v)) = (1+s_2)(1+s_1)(1+s_1)v = (1+s_2)(2)(1+s_1)v  [since (1+s_1)^2 = 2(1+s_1)]

  Actually: Φ_j = ∏_{k=1}^{j} block_k = ∏_{k=1}^{j} ∏_{m=1}^{k} (1+s_m)
           = (1+s_1)^j (1+s_2)^{j-1} ... (1+s_j)^1

  Hmm, but each block_k contains (1+s_m) for m=1,...,k.
  So in Φ_j = block_1 · block_2 · ... · block_j:
    (1+s_1) appears in ALL j blocks: factor (1+s_1)^j
    (1+s_2) appears in blocks 2,...,j: factor (1+s_2)^{j-1}
    ...
    (1+s_j) appears in block j only: factor (1+s_j)^1

  Since (1+s_m)^2 = 2(1+s_m) and (1+s_m)^k = 2^{k-1}(1+s_m) for k≥1:
    Φ_j = 2^{j-1}(1+s_1) · 2^{j-2}(1+s_2) · ... · 2^0(1+s_j)
         = 2^{j(j-1)/2} · (1+s_1)(1+s_2)...(1+s_j)

  So up to scalar, the image of Φ_j on V^H equals the image of (1+s_1)(1+s_2)...(1+s_j).

  The image after cascade j = span of {(1+s_j)(1+s_{j-1})...(1+s_1) v : v ∈ V^H}
  (the rightmost factor (1+s_1) acts first).

REVISED INTERPRETATION:
The cascade applies blocks sequentially:
After block 1: image = (1+s_1)(V^H)
After block 2: image = block_2(image after block 1) = (1+s_2)(1+s_1) applied to image after block 1
             = (1+s_2)(1+s_1)(1+s_1)(V^H) = 2(1+s_2)(1+s_1)(V^H)  [since (1+s_1)^2=2(1+s_1)]
...
So after j blocks: image ∝ (1+s_j)...(1+s_2)(1+s_1)(V^H)

This is what the code does: iteratively apply (1+s_j) to current basis.
The check is: does s_j fix the image after applying (1+s_j)?

KEY OBSERVATION from v2/v3 runs:
- For n=4,5,6,7: ALL j (odd and even) give YES: s_j fixes the image.
- dim(image) = dim(V^H) throughout (no collapsing!).
- image=V^H only at j=1 (first block); after that image ≠ V^H but s_j still fixes it.

Let's also check:
- Does s_{j+1} also fix the image? (would be interesting)
- Does EVERY s_k fix the image? (would mean image ⊆ ∩_k +1eig(s_k))
"""

from fractions import Fraction
from itertools import permutations
import time
import sys
import math

def swap_pos(sigma, i):
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
    """Gaussian elimination. Returns reduced row echelon basis."""
    rows = []
    pivot_at = {}

    for v in vecs:
        if not v:
            continue
        current = dict(v)

        # Reduce against pivots
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

def apply_block_op(basis, i):
    """Apply (I + s_i) to each vector in basis, return raw (unreduced) result."""
    result = []
    for v in basis:
        sv = apply_si_to_vec(v, i)
        new_v = vec_add(v, sv)
        if not is_zero(new_v):
            result.append(new_v)
    return result

def check_si_fixes(basis, i):
    """Check s_i v = v for all basis vectors v."""
    for v in basis:
        sv = apply_si_to_vec(v, i)
        diff = vec_sub(sv, v)
        if not is_zero(diff):
            return False
    return True

def compute_VH(n, perms_sorted, ordering):
    """V^H = joint +1 eigenspace of s_0, s_2, s_4, ... (0-indexed)."""
    H_gens = list(range(0, n-1, 2))
    current = [{sigma: Fraction(1)} for sigma in perms_sorted]
    for hi in H_gens:
        raw = apply_block_op(current, hi)
        current = reduce_basis(raw, ordering)
    return current

def run_test(n):
    print(f"\n{'='*65}")
    print(f"n = {n},  |S_n| = {math.factorial(n)}")

    t0 = time.time()
    perms_sorted = sorted(permutations(range(n)))
    ordering = perms_sorted

    VH = compute_VH(n, perms_sorted, ordering)
    dim_VH = len(VH)
    expected = math.factorial(n) // (2 ** (n // 2))
    print(f"dim(V^H) = {dim_VH}  (expected n!/2^⌊n/2⌋ = {expected}, match={dim_VH==expected})")

    if dim_VH == 0:
        print("V^H trivial, skip.")
        return

    # --- Check which generators fix V^H itself ---
    print("\nWhich s_k fix V^H itself?")
    for k in range(1, n):  # s_1,...,s_{n-1} (1-indexed) = 0,...,n-2 (0-indexed)
        fixed = check_si_fixes(VH, k-1)
        print(f"  s_{k} fixes V^H: {'YES' if fixed else 'NO'}", end="")
    print()

    print(f"\nStaircase cascade (applying (I+s_j) sequentially):")
    print(f"{'j':>4} {'par':>5} {'dim':>5} {'s_j fix':>9}", end="")
    for k in range(1, min(n, 6)):  # check s_1,...,s_5 additionally
        print(f" {'s_'+str(k):>7}", end="")
    print()
    print("-" * (30 + 8*min(n-1, 5)))

    current_basis = [dict(v) for v in VH]

    for j in range(1, n):  # j = 1,...,n-1 (1-indexed)
        i = j - 1

        raw = apply_block_op(current_basis, i)
        current_basis = reduce_basis(raw, ordering)
        dim_image = len(current_basis)
        parity = "odd" if j % 2 == 1 else "even"

        fixed_j = check_si_fixes(current_basis, i)

        print(f"  j={j}  {parity:>5}  {dim_image:5}  {'YES' if fixed_j else 'NO':>9}", end="")

        # Check which OTHER s_k fix the image
        for k in range(1, min(n, 6)):
            if k - 1 == i:
                print(f" {'(*)':>7}", end="")  # already checked
            else:
                fk = check_si_fixes(current_basis, k-1)
                print(f" {'YES' if fk else 'no':>7}", end="")
        print()
        sys.stdout.flush()

    t1 = time.time()
    print(f"\n[n={n} in {t1-t0:.2f}s]")

if __name__ == "__main__":
    print("Staircase cascade eigenspace test - FINAL")
    print("All computations use exact Fraction arithmetic")
    print()
    print("CLAIM: Φ_j(V^H) ⊆ +1 eigenspace of s_j, for all j")
    print("Checking BOTH odd and even j.")

    for n in [4, 5, 6]:
        run_test(n)
        sys.stdout.flush()

    print("\n--- n=7 (may take ~2 min) ---")
    run_test(7)
