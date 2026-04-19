"""
Staircase cascade eigenspace test.
Tests: after applying Φ_j = block_1 · ... · block_j to V^H,
does s_j fix every vector in the image?

Works in the regular representation using exact Fraction arithmetic.
s_i acts by LEFT multiplication: s_i · σ = (s_i ∘ σ).
"""

from fractions import Fraction
from itertools import permutations
import sys

def compose(sigma, tau):
    """sigma ∘ tau: apply tau first, then sigma. Both are 0-indexed tuples."""
    return tuple(sigma[tau[i]] for i in range(len(sigma)))

def swap_adjacent(sigma, i):
    """Apply s_i (swap positions i and i+1, 0-indexed) on LEFT: returns s_i ∘ sigma."""
    lst = list(sigma)
    lst[i], lst[i+1] = lst[i+1], lst[i]
    return tuple(lst)

def build_perm_list(n):
    """Return sorted list of all permutations of (0,...,n-1) and index map."""
    perms = sorted(permutations(range(n)))
    idx = {p: i for i, p in enumerate(perms)}
    return perms, idx

def build_si_matrix(perms, idx, i):
    """Build matrix of s_i in regular rep (left multiplication). Returns list of rows."""
    N = len(perms)
    M = [[Fraction(0)] * N for _ in range(N)]
    for j, sigma in enumerate(perms):
        target = swap_adjacent(sigma, i)
        k = idx[target]
        M[k][j] = Fraction(1)
    return M

def mat_vec(M, v):
    """Matrix times vector, exact."""
    n = len(M)
    return [sum(M[r][c] * v[c] for c in range(n)) for r in range(n)]

def mat_add(A, B):
    n = len(A)
    return [[A[r][c] + B[r][c] for c in range(n)] for r in range(n)]

def mat_scale(A, s):
    n = len(A)
    return [[A[r][c] * s for c in range(n)] for r in range(n)]

def identity(n):
    return [[Fraction(1) if r == c else Fraction(0) for c in range(n)] for r in range(n)]

def mat_mat(A, B):
    """Matrix multiplication, exact."""
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[Fraction(0)] * m for _ in range(n)]
    for r in range(n):
        for j in range(k):
            if A[r][j] == 0:
                continue
            for c in range(m):
                C[r][c] += A[r][j] * B[j][c]
    return C

def null_space_exact(M):
    """
    Compute null space of M over Q using exact fraction arithmetic.
    Returns list of basis vectors for ker(M).
    M is a list of rows (list of lists of Fraction).
    """
    rows = [row[:] for row in M]
    nrows = len(rows)
    ncols = len(rows[0]) if nrows > 0 else 0

    pivot_cols = []
    pivot_row = 0

    for col in range(ncols):
        # Find pivot
        found = -1
        for row in range(pivot_row, nrows):
            if rows[row][col] != 0:
                found = row
                break
        if found == -1:
            continue
        rows[pivot_row], rows[found] = rows[found], rows[pivot_row]
        pivot = rows[pivot_row][col]
        rows[pivot_row] = [x / pivot for x in rows[pivot_row]]
        for row in range(nrows):
            if row != pivot_row and rows[row][col] != 0:
                factor = rows[row][col]
                rows[row] = [rows[row][c] - factor * rows[pivot_row][c] for c in range(ncols)]
        pivot_cols.append(col)
        pivot_row += 1

    free_cols = [c for c in range(ncols) if c not in pivot_cols]

    basis = []
    for fc in free_cols:
        vec = [Fraction(0)] * ncols
        vec[fc] = Fraction(1)
        for i, pc in enumerate(pivot_cols):
            if pivot_row > i:
                vec[pc] = -rows[i][fc]
        basis.append(vec)

    return basis

def eigenspace_plus1(M_si, N):
    """
    Compute +1 eigenspace of matrix M_si (N x N).
    = null space of (M_si - I).
    """
    diff = [[M_si[r][c] - (Fraction(1) if r == c else Fraction(0)) for c in range(N)] for r in range(N)]
    return null_space_exact(diff)

def intersect_eigenspaces(bases, N):
    """
    Intersect +1 eigenspaces given by their bases.
    Each basis is a list of vectors in Q^N.
    Intersection = vectors in span of first basis that also lie in span of second, etc.
    We do this by sequential intersection.
    """
    if not bases:
        return [standard_basis_vec(i, N) for i in range(N)]

    current = bases[0]
    for b in bases[1:]:
        current = intersect_two_spaces(current, b, N)
    return current

def intersect_two_spaces(basis1, basis2, N):
    """
    Intersect two subspaces given by their bases.
    Returns basis for intersection.
    """
    if not basis1 or not basis2:
        return []

    k1 = len(basis1)
    k2 = len(basis2)

    # A vector in intersection: x = sum a_i * b1_i = sum c_j * b2_j
    # => sum a_i * b1_i - sum c_j * b2_j = 0
    # Build matrix [B1 | -B2] and find null space
    # Columns are basis vectors

    mat = []
    for row in range(N):
        r = [basis1[i][row] for i in range(k1)] + [-basis2[j][row] for j in range(k2)]
        mat.append(r)

    null = null_space_exact(mat)

    # Each null vector (a1,...,ak1, c1,...,ck2) gives intersection vector sum a_i * b1_i
    result = []
    for nv in null:
        coeffs = nv[:k1]
        vec = [sum(coeffs[i] * basis1[i][row] for i in range(k1)) for row in range(N)]
        # Verify it's nonzero
        if any(v != 0 for v in vec):
            result.append(vec)

    return result

def apply_matrix_to_basis(M, basis):
    """Apply matrix M to each basis vector."""
    return [mat_vec(M, v) for v in basis]

def vectors_in_span(vecs, basis, N):
    """
    Check if each vector in vecs lies in span of basis.
    Returns True/False.
    """
    if not basis:
        return all(all(v == 0 for v in vec) for vec in vecs)
    if not vecs:
        return True

    k = len(basis)

    for vec in vecs:
        # Solve: sum c_i * basis[i] = vec
        # Augmented matrix [B | vec]
        mat = []
        for row in range(N):
            r = [basis[i][row] for i in range(k)] + [vec[row]]
            mat.append(r)

        # Row reduce
        rows = [r[:] for r in mat]
        nrows = len(rows)
        ncols = len(rows[0])

        pivot_row = 0
        consistent = True
        for col in range(ncols - 1):
            found = -1
            for row in range(pivot_row, nrows):
                if rows[row][col] != 0:
                    found = row
                    break
            if found == -1:
                continue
            rows[pivot_row], rows[found] = rows[found], rows[pivot_row]
            pivot = rows[pivot_row][col]
            rows[pivot_row] = [x / pivot for x in rows[pivot_row]]
            for row in range(nrows):
                if row != pivot_row and rows[row][col] != 0:
                    factor = rows[row][col]
                    rows[row] = [rows[row][c] - factor * rows[pivot_row][c] for c in range(ncols)]
            pivot_row += 1

        # Check for inconsistency: row with all zeros in coefficient part but nonzero in rhs
        for row in range(nrows):
            if all(rows[row][c] == 0 for c in range(ncols - 1)) and rows[row][-1] != 0:
                consistent = False
                break

        if not consistent:
            return False

    return True

def is_fixed_by(M_si, image_basis, N):
    """
    Check if s_i fixes every vector in the image (given by basis).
    Equivalent: M_si * v = v for all v in image.
    Equivalent: (M_si - I) * v = 0 for all v in image.
    Equivalent: image_basis ⊆ +1 eigenspace of M_si.
    """
    if not image_basis:
        return True  # Empty image is vacuously fixed

    for v in image_basis:
        Mv = mat_vec(M_si, v)
        diff = [Mv[i] - v[i] for i in range(N)]
        if any(d != 0 for d in diff):
            return False
    return True

def run_test(n):
    print(f"\n{'='*60}")
    print(f"n = {n}, |S_n| = {1}")

    import math
    fact_n = math.factorial(n)
    print(f"n = {n}, |S_n| = {fact_n}")

    perms, idx = build_perm_list(n)
    N = len(perms)

    # Build s_i matrices for i = 0, ..., n-2 (0-indexed transpositions)
    print(f"Building {n-1} transposition matrices ({N}x{N})...")
    Si = []
    for i in range(n-1):
        Si.append(build_si_matrix(perms, idx, i))

    # Compute V^H = joint +1 eigenspace of s_0, s_2, s_4, ... (0-indexed: generators of H)
    # H = <s_1, s_3, s_5, ...> in 1-indexed = <s_0, s_2, s_4, ...> in 0-indexed
    H_gens = list(range(0, n-1, 2))  # 0-indexed: s_0, s_2, s_4, ...
    print(f"H generators (0-indexed): {H_gens} = s_{{{', '.join(str(i+1) for i in H_gens)}}} in 1-indexed")

    # Compute +1 eigenspaces for each H generator and intersect
    h_bases = []
    for i in H_gens:
        basis_i = eigenspace_plus1(Si[i], N)
        h_bases.append(basis_i)
        print(f"  dim(+1 eigenspace of s_{i+1}) = {len(basis_i)}")

    if len(h_bases) == 1:
        VH = h_bases[0]
    elif len(h_bases) > 1:
        VH = h_bases[0]
        for b in h_bases[1:]:
            VH = intersect_two_spaces(VH, b, N)
    else:
        VH = [([Fraction(1) if j == k else Fraction(0) for j in range(N)]) for k in range(N)]

    print(f"dim(V^H) = {len(VH)}")

    if len(VH) == 0:
        print("V^H is trivial, skipping.")
        return

    # Now apply staircase blocks one by one
    # block_j (1-indexed j) = (I + S_{j-1}) (0-indexed: S_{j-1})
    # Φ_j = block_1 · block_2 · ... · block_j
    # Applied to v: first apply block_1, then block_2, etc.
    # So image after j blocks = block_j applied to (image after j-1 blocks)

    # P_i = (I + Si[i]) / 2  is the projector, but we use (I + Si[i]) (unnormalized)
    # The claim is about the IMAGE, so scaling doesn't matter for eigenspace membership.

    print(f"\nStaircase cascade results:")
    print(f"{'Block j':>8} {'j-parity':>10} {'dim(image)':>12} {'s_j fixes image?':>18}")
    print("-" * 55)

    current_basis = VH[:]

    for j in range(1, n):  # j = 1, ..., n-1 (1-indexed blocks)
        # block_j = (I + S_{j-1}) (0-indexed transposition index = j-1)
        i = j - 1  # 0-indexed

        # Apply (I + Si[i]) to each vector in current_basis
        new_basis_raw = []
        for v in current_basis:
            Siv = mat_vec(Si[i], v)
            new_v = [v[k] + Siv[k] for k in range(N)]
            new_basis_raw.append(new_v)

        # The image of V^H under Φ_j is the span of new_basis_raw
        # But we need to find a basis (reduce to linearly independent set)
        current_basis = reduce_to_basis(new_basis_raw, N)

        dim_image = len(current_basis)

        # Check: does s_j fix every vector in the image?
        # s_j in 1-indexed = Si[j-1] in 0-indexed = Si[i]
        fixed = is_fixed_by(Si[i], current_basis, N)

        parity = "odd" if j % 2 == 1 else "even"
        fixed_str = "YES" if fixed else "NO"

        print(f"  j={j:2d}   {parity:>8}   dim={dim_image:4d}   s_{j} fixes image: {fixed_str}")

def reduce_to_basis(vecs, N):
    """
    Given a list of vectors in Q^N, return a basis for their span.
    Uses Gaussian elimination.
    """
    if not vecs:
        return []

    # Build matrix with vecs as rows
    rows = [v[:] for v in vecs]
    nrows = len(rows)
    ncols = N

    pivot_cols = []
    pivot_row = 0

    for col in range(ncols):
        found = -1
        for row in range(pivot_row, nrows):
            if rows[row][col] != 0:
                found = row
                break
        if found == -1:
            continue
        rows[pivot_row], rows[found] = rows[found], rows[pivot_row]
        pivot = rows[pivot_row][col]
        rows[pivot_row] = [x / pivot for x in rows[pivot_row]]
        for row in range(nrows):
            if row != pivot_row and rows[row][col] != 0:
                factor = rows[row][col]
                rows[row] = [rows[row][c] - factor * rows[pivot_row][c] for c in range(ncols)]
        pivot_cols.append(col)
        pivot_row += 1

    # Return only the pivot rows
    return [rows[i] for i in range(len(pivot_cols))]

if __name__ == "__main__":
    import time

    print("Staircase cascade eigenspace test")
    print("Claim: After Φ_j applied to V^H, image lies in +1 eigenspace of s_j")
    print("Using EXACT Fraction arithmetic in regular representation")

    for n in [4, 5, 6]:
        t0 = time.time()
        run_test(n)
        t1 = time.time()
        print(f"[n={n} took {t1-t0:.1f}s]")

    # n=7 might be slow (5040 x 5040), try if time permits
    print("\nAttempting n=7 (5040x5040 matrices)...")
    print("This may take a while with exact arithmetic...")
    t0 = time.time()
    run_test(7)
    t1 = time.time()
    print(f"[n=7 took {t1-t0:.1f}s]")
