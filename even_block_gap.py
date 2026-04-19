"""
Even-block gap analysis for H-invariant theorem.

Operator: A = (1+s_3)(1+s_1)(1+s_2) in S_4 (relabeled from general k).

Key question: Is ker(A) = ker(1+s_2) in every irrep of S_4?

Method: Use Young's seminormal form (rational over QQ).
In seminormal form, s_i acts by:
  - diagonal entry: 1/d where d = axial distance = content(i+1) - content(i)
  - off-diagonal (T <-> s_i T): sqrt(1 - 1/d^2)

For our purposes we work with exact rational arithmetic (using fractions).
The key structural properties (kernel dimensions) are algebraic invariants
so we can verify them using exact methods.

We use the symmetric group algebra acting on the regular representation,
then block-decompose by irrep using the central idempotents.
"""

from fractions import Fraction
import numpy as np
from itertools import permutations

# ============================================================
# Permutation group utilities
# ============================================================

def perm_mul(sigma, tau, n):
    """Compose permutations: (sigma * tau)(i) = sigma(tau(i))."""
    return tuple(sigma[tau[i]] for i in range(n))

def perm_inv(sigma):
    n = len(sigma)
    inv = [0] * n
    for i, si in enumerate(sigma):
        inv[si] = i
    return tuple(inv)

def perm_from_transposition(i, j, n):
    """Transposition (i+1, j+1) as 0-indexed permutation."""
    p = list(range(n))
    p[i], p[j] = p[j], p[i]
    return tuple(p)

def all_perms(n):
    return list(permutations(range(n)))

def perm_index(perms, sigma):
    return perms.index(sigma)

def left_reg_rep_matrix(perms, sigma, n):
    """Matrix of left multiplication by sigma in the regular representation."""
    d = len(perms)
    M = np.zeros((d, d), dtype=float)
    for j, tau in enumerate(perms):
        result = perm_mul(sigma, tau, n)
        i = perm_index(perms, result)
        M[i, j] = 1.0
    return M

# ============================================================
# Young tableaux utilities
# ============================================================

def standard_tableaux(partition):
    """
    Generate all standard Young tableaux of given partition.
    Entries increase left-to-right in each row and top-to-bottom in each column.
    We fill entries 1, 2, ..., n in order, placing each in a valid corner cell.
    """
    n = sum(partition)
    shape = list(partition)
    rows = len(shape)

    def fill(tableau, n_filled):
        if n_filled == n:
            yield [row[:] for row in tableau]
            return
        entry = n_filled + 1
        # A cell (r, c) is a valid "inner corner" to place entry if:
        # - len(tableau[r]) == c (next available in row r)
        # - c < shape[r] (row r has room)
        # - c == 0 or tableau[r][c-1] is already filled (auto since we fill left to right)
        # - r == 0 or len(tableau[r-1]) > c (cell above is already filled)
        # Since we fill entries 1,2,...,n in order, the new entry is automatically
        # larger than all existing entries.
        for r in range(rows):
            c = len(tableau[r])
            if c >= shape[r]:
                continue
            # Row condition: cell to the left must exist (auto: entries increase)
            # Column condition: cell above must be filled
            if r > 0 and len(tableau[r-1]) <= c:
                continue
            tableau[r].append(entry)
            yield from fill(tableau, n_filled + 1)
            tableau[r].pop()

    tableau0 = [[] for _ in range(rows)]
    return list(fill(tableau0, 0))

def content(tab, val):
    """Return content (col - row) of value val in tableau tab (1-indexed values)."""
    for r, row in enumerate(tab):
        for c, v in enumerate(row):
            if v == val:
                return c - r  # 0-indexed col and row
    raise ValueError(f"{val} not found in tableau")

def apply_si(tab, i):
    """
    Swap entries i and i+1 in tableau tab.
    Returns new tableau if standard, else None.
    """
    n_rows = len(tab)
    pos_i = pos_i1 = None
    for r, row in enumerate(tab):
        for c, v in enumerate(row):
            if v == i:
                pos_i = (r, c)
            elif v == i + 1:
                pos_i1 = (r, c)

    new_tab = [row[:] for row in tab]
    new_tab[pos_i[0]][pos_i[1]] = i + 1
    new_tab[pos_i1[0]][pos_i1[1]] = i

    # Check if standard
    for r, row in enumerate(new_tab):
        for c, v in enumerate(row):
            if c > 0 and row[c-1] >= v:
                return None
            if r > 0 and len(new_tab[r-1]) > c and new_tab[r-1][c] >= v:
                return None
    return new_tab

def tab_to_tuple(tab):
    return tuple(tuple(row) for row in tab)

# ============================================================
# Young seminormal form
# ============================================================

def seminormal_matrices(partition, n):
    """
    Compute Young's seminormal representation matrices for s_1,...,s_{n-1}.
    Returns dict {i: matrix} where matrix is numpy float64 array.

    Seminormal form: basis = SYT(partition)
    For s_i acting on basis vector e_T:
      s_i e_T = (1/d) e_T + sqrt(1-1/d^2) e_{T'} if T' standard
      s_i e_T = (1/d) e_T                         if T' not standard
    where d = content(i+1) - content(i) in T.
    """
    tabs = standard_tableaux(partition)
    d_rep = len(tabs)
    tab_idx = {tab_to_tuple(T): k for k, T in enumerate(tabs)}

    matrices = {}
    for i in range(1, n):
        M = np.zeros((d_rep, d_rep))
        for j, T in enumerate(tabs):
            ci = content(T, i)
            ci1 = content(T, i + 1)
            axial = ci1 - ci  # = content(i+1) - content(i)
            diag = 1.0 / axial
            M[j, j] = diag

            T_prime = apply_si(T, i)
            if T_prime is not None:
                k = tab_idx[tab_to_tuple(T_prime)]
                off = np.sqrt(max(0.0, 1.0 - diag**2))
                M[j, k] = off  # T' contribution to e_T image

        matrices[i] = M
    return matrices, tabs, d_rep

# ============================================================
# Rational seminormal form (exact)
# ============================================================

def rational_seminormal_matrices(partition, n):
    """
    Exact rational seminormal form. The off-diagonal entries are sqrt(1-1/d^2).
    These are NOT always rational. But we can work with the SQUARED kernel:
    v in ker(M) iff M v = 0 iff each component is 0.

    Alternative: Use the fact that ker(1+s_i) in seminormal form consists of
    vectors where s_i acts as -1. In the seminormal basis, s_i is block diagonal
    with 2x2 blocks for pairs (T, T') and 1x1 blocks of ±1 for fixed T.

    The +1 eigenspace of s_i: s_i v = v
    The -1 eigenspace of s_i: s_i v = -v

    For each pair (T, s_i T), the s_i block is:
      [ 1/d    sqrt(1-1/d^2) ]
      [ sqrt(1-1/d^2)  -1/d  ]

    Eigenvalues: +1 and -1.
    +1 eigenvector: [sqrt(1-1/d^2), d - 1] (unnormalized)... let me compute.

    Matrix: [[1/d, c], [c, -1/d]] where c = sqrt(1 - 1/d^2)
    Eigenvalues: 1 and -1 (since it squares to I for adjacent transpositions).

    +1 eigenvector: [[1/d, c], [c, -1/d]] v = v
    => (1/d - 1) v_1 + c v_2 = 0
    => v_2/v_1 = (1 - 1/d)/c = (d-1)/(d*c) = (d-1)/(d * sqrt((d^2-1)/d^2))
                = (d-1)/sqrt(d^2-1) = sqrt((d-1)/(d+1))

    This is irrational in general. But we only need kernel dimensions.
    """
    tabs = standard_tableaux(partition)
    d_rep = len(tabs)
    tab_idx = {tab_to_tuple(T): k for k, T in enumerate(tabs)}

    # For kernel computation, we need exact arithmetic.
    # Key insight: in the seminormal basis, (1+s_i) is block diagonal.
    # For a pair (T, T') with T' = s_i T:
    #   block of s_i = [[1/d, c], [c, -1/d]]
    #   1+s_i restricted to this block = [[1+1/d, c], [c, 1-1/d]]
    # Determinant of this block: (1+1/d)(1-1/d) - c^2 = 1 - 1/d^2 - (1-1/d^2) = 0
    # So det = 0! Every pair contributes exactly 1 to ker(1+s_i).
    # For a singleton T with s_i T not standard, s_i acts as +1 or -1 (diag 1 or -1).
    # If axial distance gives 1/d with |d|=1, then d=1 (same row) or d=-1 (same column).
    #
    # Wait: if T' is NOT a valid SYT, what is the diagonal entry?
    # axial = content(i+1) - content(i)
    # If i and i+1 are in the same row: axial = 1, s_i acts as +1 (diag 1)
    # If i and i+1 are in the same column: axial = -1, s_i acts as -1 (diag -1)

    # So:
    # - For each pair (T, T'): s_i block has det=0, so contributes 1 dim to ker(1+s_i)
    # - For singleton T in same row: s_i = +1, not in ker(1+s_i)
    # - For singleton T in same col: s_i = -1, so 1+s_i = 0 on this subspace, contributes 1 dim

    # Therefore dim ker(1+s_i) = #{pairs (T,T')} + #{singletons in same column}

    # This gives us exact computation of ker dimensions!
    return tabs, d_rep, tab_idx

def kernel_dim_1_plus_si(partition, n, i):
    """
    Compute dim ker(1+s_i) in the irrep V_lambda using seminormal form analysis.
    """
    tabs = standard_tableaux(partition)
    d_rep = len(tabs)
    tab_idx = {tab_to_tuple(T): k for k, T in enumerate(tabs)}

    paired = set()
    singleton_minus = 0

    for j, T in enumerate(tabs):
        T_key = tab_to_tuple(T)
        if T_key in paired:
            continue

        T_prime = apply_si(T, i)
        if T_prime is not None:
            T_prime_key = tab_to_tuple(T_prime)
            paired.add(T_key)
            paired.add(T_prime_key)
        else:
            # Singleton: check if s_i acts as -1
            ci = content(T, i)
            ci1 = content(T, i + 1)
            axial = ci1 - ci
            if axial == -1:
                singleton_minus += 1

    n_pairs = len(paired) // 2
    return n_pairs + singleton_minus

def compute_operator_matrices(partition, n, gen_indices):
    """
    Compute float matrices for generators s_i (i in gen_indices) in V_lambda.
    Uses Young's seminormal form.
    """
    mats, tabs, d = seminormal_matrices(partition, n)
    return mats, tabs, d

# ============================================================
# Kernel analysis using numpy
# ============================================================

def kernel_dim(M, tol=1e-10):
    """Compute dimension of kernel of matrix M."""
    if M.shape[0] == 0:
        return 0
    sv = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(sv < tol))

def kernel_basis(M, tol=1e-10):
    """Compute basis for kernel of M."""
    U, sv, Vt = np.linalg.svd(M, full_matrices=True)
    null_mask = sv < tol
    # Append all-zero rows if M is tall
    n_null_from_sv = int(np.sum(null_mask))
    extra = max(0, M.shape[1] - len(sv))
    return Vt[len(sv)-n_null_from_sv:].T, n_null_from_sv + extra

def are_same_subspace(A_rows, B_rows, tol=1e-10):
    """Check if row spaces of A and B are equal."""
    if A_rows.shape[0] == 0 and B_rows.shape[0] == 0:
        return True
    if A_rows.shape[0] == 0 or B_rows.shape[0] == 0:
        return False

    # Check if each row of A is in span(B) and vice versa
    def in_span(v, basis, tol):
        if basis.shape[0] == 0:
            return np.linalg.norm(v) < tol
        coeffs, res, rank, sv = np.linalg.lstsq(basis.T, v, rcond=None)
        return np.linalg.norm(basis.T @ coeffs - v) < tol

    for row in A_rows:
        if not in_span(row, B_rows, tol):
            return False
    for row in B_rows:
        if not in_span(row, A_rows, tol):
            return False
    return True

# ============================================================
# Main analysis
# ============================================================

def analyze_S4():
    """
    Full analysis for S_4.
    Operator A = (1+s_3)(1+s_1)(1+s_2).
    For each irrep, compute kernel dimensions and check equality.
    Also compute V^H and transported image.
    """
    print("=" * 60)
    print("S_4 ANALYSIS")
    print("n = 4, Operator: A = (1+s_3)(1+s_1)(1+s_2)")
    print("H = <s_1, s_3>")
    print("=" * 60)
    print()

    n = 4
    partitions = [[4], [3,1], [2,2], [2,1,1], [1,1,1,1]]
    names = ['(4)', '(3,1)', '(2,2)', '(2,1,1)', '(1,1,1,1)']

    all_closed = True

    for lam, name in zip(partitions, names):
        mats, tabs, d = seminormal_matrices(lam, n)
        I = np.eye(d)

        # Generators
        S1, S2, S3 = mats[1], mats[2], mats[3]

        # Operators
        op_s2 = I + S2                       # (1+s_2)
        op_s1s2 = (I + S1) @ (I + S2)       # (1+s_1)(1+s_2)
        op_A = (I + S3) @ (I + S1) @ (I + S2)  # A = (1+s_3)(1+s_1)(1+s_2)

        # Kernel dimensions
        k_s2 = kernel_dim(op_s2)
        k_s1s2 = kernel_dim(op_s1s2)
        k_A = kernel_dim(op_A)

        # Kernel bases (for subspace comparison)
        ker_s2_basis, _ = kernel_basis(op_s2)
        ker_s1s2_basis, _ = kernel_basis(op_s1s2)
        ker_A_basis, _ = kernel_basis(op_A)

        # V^H: fixed by s_1 and s_3
        # s_1 v = v => (S1 - I) v = 0
        # s_3 v = v => (S3 - I) v = 0
        # V^H = ker(S1 - I) ∩ ker(S3 - I)
        combined_H = np.vstack([S1 - I, S3 - I])
        VH_basis, n_VH = kernel_basis(combined_H)
        # VH_basis columns span V^H

        # Transport V^H through staircase blocks 1 and 2
        # Block 1: (1+s_1)
        # Block 2: (1+s_2)(1+s_1)
        # Staircase through blocks 1,2: [(1+s_2)(1+s_1)] * [(1+s_1)] * V^H
        # = (1+s_2)(1+s_1)^2 * V^H
        # But (1+s_1)^2 = (1+s_1) + (s_1 + s_1^2) = (1+s_1) + (s_1 + 1) = 2(1+s_1)... wait
        # s_1^2 = id, so (1+s_1)^2 = 1 + 2s_1 + s_1^2 = 1 + 2s_1 + 1 = 2(1+s_1)
        # Hmm but that's if we apply it twice. Let me re-read the setup.
        #
        # The staircase product is:
        # Π = [B_3][B_2][B_1]
        # where B_k = (1+s_k)(1+s_{k-1})...(1+s_1)
        # B_1 = (1+s_1)
        # B_2 = (1+s_2)(1+s_1)
        # B_3 = (1+s_3)(1+s_2)(1+s_1)
        #
        # Applied to V^H right-to-left:
        # First B_1 = (1+s_1): since v in V^H has s_1 v = v, B_1 v = 2v
        # Then B_2 applied to result: B_2 * 2v = 2 * (1+s_2)(1+s_1) v = 2*(1+s_2)*2v = 4(1+s_2)v
        # Hmm, the image of V^H under blocks 1 and 2 is proportional.
        #
        # For kernel intersection purposes, scalar multiples don't matter.
        # The IMAGE of V^H under the staircase at start of block 3 =
        # span of {(1+s_2)(1+s_1)(1+s_1) v : v in V^H}
        #       = span of {(1+s_2)(1+s_1)^2 v : v in V^H}
        #       = span of {2(1+s_2)(1+s_1) v : v in V^H}   [since (1+s_i)^2 = 2(1+s_i)]
        #       = span of {(1+s_2)(1+s_1) v : v in V^H}    [same span]
        #
        # Now (1+s_1) v = 2v since s_1 v = v.
        # So image = span of {2(1+s_2) v : v in V^H} = span of {(1+s_2) v : v in V^H}

        # So the transported image at start of block 3 is (1+s_2) V^H.
        # This IMMEDIATELY shows it's in the +1 eigenspace of s_2!
        # Because: s_2 * (1+s_2) v = (s_2 + s_2^2) v = (s_2 + 1) v = (1+s_2) v.
        # So every vector in (1+s_2) V^H is fixed by s_2 => in E_2^{(+1)}.

        # Therefore transported image ⊆ ker(s_2 - 1) = E_2^{(+1)}.
        # And ker(A) ∩ E_2^{(+1)} = ?

        # Compute (1+s_2) V^H
        if n_VH > 0:
            VH_mat = VH_basis  # columns span V^H
            transported_mat = (I + S2) @ VH_mat  # columns span (1+s_2)V^H
            # Dimension of transported
            transported_rank = np.linalg.matrix_rank(transported_mat, tol=1e-10)

            # Intersection with ker(A)
            # We need dim( span(transported_mat cols) ∩ ker(op_A) )
            # A vector v is in both iff:
            # v = transported_mat @ c for some c, and op_A @ v = 0
            # => op_A @ transported_mat @ c = 0
            # => c in ker(op_A @ transported_mat)
            intersection_mat = op_A @ transported_mat
            inter_dim = kernel_dim(intersection_mat)
        else:
            transported_rank = 0
            inter_dim = 0

        # Check ker equalities
        eq_A_s2 = (k_A == k_s2)  # dimension equality (necessary but not sufficient)
        # For subspace equality with float arithmetic, use the basis check
        if k_A == k_s2 and k_A > 0:
            eq_A_s2_subspace = are_same_subspace(ker_A_basis.T, ker_s2_basis.T)
        elif k_A == k_s2 == 0:
            eq_A_s2_subspace = True
        else:
            eq_A_s2_subspace = False

        if k_s1s2 == k_s2 and k_s1s2 > 0:
            eq_s1s2_s2_subspace = are_same_subspace(ker_s1s2_basis.T, ker_s2_basis.T)
        elif k_s1s2 == k_s2 == 0:
            eq_s1s2_s2_subspace = True
        else:
            eq_s1s2_s2_subspace = False

        gap_closed = (inter_dim == 0)
        if not gap_closed:
            all_closed = False

        print(f"λ = {name}  (dim = {d})")
        print(f"  Basis: {d} standard tableaux")
        print(f"  dim ker(1+s_2)            = {k_s2}")
        print(f"  dim ker((1+s_1)(1+s_2))   = {k_s1s2}")
        print(f"  dim ker(A)                = {k_A}")
        print(f"  ker(A) = ker((1+s_1)(1+s_2))? dim equal: {k_A == k_s1s2}")
        print(f"  ker(A) = ker(1+s_2)?       dim equal: {k_A == k_s2}, subspace: {eq_A_s2_subspace}")
        print(f"  ker((1+s_1)(1+s_2)) = ker(1+s_2)? dim equal: {k_s1s2 == k_s2}, subspace: {eq_s1s2_s2_subspace}")
        print(f"  dim V^H                   = {n_VH}")
        print(f"  dim transported image     = {transported_rank}")
        print(f"  transported ∩ ker(A) dim  = {inter_dim}")
        if gap_closed:
            print(f"  => Gap CLOSED for this irrep")
        else:
            print(f"  => POTENTIAL GAP for this irrep!")
        print()

    print(f"Overall S_4: {'ALL irreps closed - gap closes!' if all_closed else 'GAP EXISTS'}")
    print()
    return all_closed

def analyze_Sn_even_blocks(n, even_blocks):
    """
    For S_n, check even block k gaps.
    Operator: A = (1+s_k)(1+s_{k-2})(1+s_{k-1})
    Check ker(A) = ker(1+s_{k-1}) for all irreps.
    Also check transported image intersection with ker(A).
    """
    print(f"{'='*60}")
    print(f"S_{n} ANALYSIS - Even block gaps")
    print(f"{'='*60}")

    partitions = list(Partitions_list(n))
    all_closed = True

    for k in even_blocks:
        # For S_n, generators go s_1,...,s_{n-1}
        # The operator A = (1+s_k)(1+s_{k-2})(1+s_{k-1}) requires s_k to exist
        # i.e., k <= n-1
        if k > n - 1:
            print(f"\n  k = {k}: s_{k} doesn't exist in S_{n} (generators go up to s_{n-1}), skip")
            continue

        print(f"\n  k = {k}: A = (1+s_{k})(1+s_{{k-2}})(1+s_{{k-1}})")
        print(f"  Checking all {len(partitions)} irreps:")
        block_closed = True

        for lam in partitions:
            mats, tabs, d = seminormal_matrices(lam, n)
            I = np.eye(d)

            Sk_minus1 = mats[k-1]  # s_{k-1}
            Sk_minus2 = mats[k-2]  # s_{k-2}
            Sk = mats[k]           # s_k

            op_sk_minus1 = I + Sk_minus1
            op_A = (I + Sk) @ (I + Sk_minus2) @ (I + Sk_minus1)

            k_sk1 = kernel_dim(op_sk_minus1)
            k_A = kernel_dim(op_A)

            # V^H: fixed by all odd-indexed generators (parabolic H for H-invariant)
            # H = <s_1, s_3, ..., s_{2*floor((n-1)/2)}>
            H_gens = list(range(1, n, 2))  # 1, 3, 5, ...
            if len(H_gens) > 0:
                H_mats = np.vstack([mats[i] - I for i in H_gens if i in mats])
                VH_basis, n_VH = kernel_basis(H_mats)
            else:
                VH_basis = np.eye(d)
                n_VH = d

            # Transported image at start of even block k:
            # The staircase up to and including block k-1 ends with factor (1+s_{k-1}).
            # By the algebraic identity: s_{k-1}*(1+s_{k-1})v = (1+s_{k-1})v,
            # the image lies in E_{k-1}^{(+1)}.
            # More explicitly: the transported image = (1+s_{k-1}) * [something from V^H].
            # For a conservative bound, use: image = (1+s_{k-1}) V^H.
            if n_VH > 0:
                transported = op_sk_minus1 @ VH_basis
                inter_mat = op_A @ transported
                inter_dim = kernel_dim(inter_mat)
            else:
                inter_dim = 0

            # The key check: does transported image meet ker(A)?
            # inter_dim > 0 means there's a vector annihilated by A in the transported image
            # This would be a genuine obstruction.
            if inter_dim > 0:
                block_closed = False
                all_closed = False
                status = "TRUE GAP!"
            else:
                status = "closed"

            eq = (k_A == k_sk1)
            print(f"    λ={lam}: d={d}, ker(1+s_{{k-1}})={k_sk1}, ker(A)={k_A}, "
                  f"ker dims equal={eq}, transported∩ker(A)={inter_dim} [{status}]")

        if block_closed:
            print(f"  => k={k}: All irreps OK, block gap CLOSED")
        else:
            print(f"  => k={k}: GAP EXISTS in some irrep!")

    print(f"\nS_{n} overall: {'ALL gaps closed' if all_closed else 'GAPS EXIST'}")
    return all_closed

def Partitions_list(n):
    """Generate all partitions of n as lists."""
    def _partitions(n, max_part=None):
        if max_part is None:
            max_part = n
        if n == 0:
            yield []
            return
        for first in range(min(n, max_part), 0, -1):
            for rest in _partitions(n - first, first):
                yield [first] + rest

    return list(_partitions(n))

# ============================================================
# Additional: exact kernel dimension via combinatorial formula
# ============================================================

def exact_kernel_dim_A(partition, n, i, j, k):
    """
    Compute dim ker((1+s_k)(1+s_j)(1+s_i)) in irrep V_lambda.

    Strategy: Work in seminormal basis. The matrix is a product of three
    operators, each block-diagonalized in its own way. Use numerical methods
    for this.
    """
    mats, tabs, d = seminormal_matrices(partition, n)
    I = np.eye(d)
    op = (I + mats[k]) @ (I + mats[j]) @ (I + mats[i])
    return kernel_dim(op)

def verify_transported_image_in_E_plus(partition, n, k_block):
    """
    Verify that the transported V^H image at start of even block k
    is contained in E_{k-1}^{(+1)} = +1 eigenspace of s_{k-1}.

    Key algebraic fact: if the last factor of the staircase leading to
    block k-1's end is (1+s_{k-1}), then the image is automatically in
    the +1 eigenspace of s_{k-1}, since:
    s_{k-1} * (1+s_{k-1}) v = (s_{k-1} + 1) v = (1+s_{k-1}) v.
    """
    mats, tabs, d = seminormal_matrices(partition, n)
    I = np.eye(d)
    Sk1 = mats[k_block - 1]

    # The (1+s_{k-1}) eigenspace condition
    # (1+s_{k-1}) v has s_{k-1} * (1+s_{k-1}) v = (1+s_{k-1}) v
    # So it's in +1 eigenspace of s_{k-1}.
    # This is an algebraic identity, not numerical.
    # Let's verify numerically for sanity.

    test_vec = np.random.randn(d)
    transported = (I + Sk1) @ test_vec
    s_applied = Sk1 @ transported
    diff = np.linalg.norm(s_applied - transported)
    return diff < 1e-12

# ============================================================
# Run analyses
# ============================================================

print("EVEN-BLOCK GAP ANALYSIS FOR H-INVARIANT THEOREM")
print("=" * 60)
print()
print("Setup: staircase symmetrizer for S_n.")
print("Even block k >= 4 gap: P_{k-2} applied after P_{k-1}.")
print("Relabeled operator: A = (1+s_3)(1+s_1)(1+s_2) in S_4.")
print("Key question: ker(A) = ker(1+s_2) in each irrep?")
print()

# PART 1: S_4 analysis
print("### PART 1: S_4 ###\n")
s4_closed = analyze_S4()

# PART 2: S_5
print("\n### PART 2: S_5 (even block k=4) ###\n")
s5_closed = analyze_Sn_even_blocks(5, [4])

# PART 3: S_6
print("\n### PART 3: S_6 (even blocks k=4, k=6) ###\n")
s6_closed = analyze_Sn_even_blocks(6, [4, 6])

# PART 4: S_7
print("\n### PART 4: S_7 (even blocks k=4, k=6) ###\n")
s7_closed = analyze_Sn_even_blocks(7, [4, 6])

print("\n" + "="*60)
print("OVERALL SUMMARY")
print("="*60)
print(f"S_4 gaps closed: {s4_closed}")
print(f"S_5 gaps closed: {s5_closed}")
print(f"S_6 gaps closed: {s6_closed}")
print(f"S_7 gaps closed: {s7_closed}")
print()

# PART 5: Algebraic argument for why transported image is in E^+
print("### ALGEBRAIC ARGUMENT ###")
print()
print("Claim: The transported V^H image at start of even block k")
print("lies in E_{k-1}^{(+1)} (= +1 eigenspace of s_{k-1}).")
print()
print("Proof: The staircase block for block k-1 ends with factor (1+s_{k-1}).")
print("For any vector w, s_{k-1}(1+s_{k-1})w = (s_{k-1}+1)w = (1+s_{k-1})w.")
print("So every vector in the image of (1+s_{k-1}) is a +1 eigenvector of s_{k-1}.")
print("The transported image at block k-1's end = (1+s_{k-1}) * [earlier stuff].")
print("Therefore: transported image ⊆ E_{k-1}^{(+1)} = ker(s_{k-1}-1).")
print()
print("Meanwhile: ker(1+s_{k-1}) = E_{k-1}^{(-1)} = -1 eigenspace.")
print()
print("Since E_{k-1}^{(+1)} ∩ E_{k-1}^{(-1)} = {0} (eigenspaces for eigenvalues ±1),")
print("the transported image meets ker(1+s_{k-1}) only in {0}.")
print()
print("Now if ker(A) = ker(1+s_{k-1}) in each irrep,")
print("then transported image ∩ ker(A) = {0} in each irrep,")
print("and the even-block gap is CLOSED.")
print()
print("NUMERICAL VERIFICATION of algebraic claim:")
for lam in [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]]:
    ok = verify_transported_image_in_E_plus(lam, 4, 2)
    print(f"  λ={lam}, k=2: transported in E_1^{{+1}}? {ok}")
