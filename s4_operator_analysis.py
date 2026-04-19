"""
S_4 Operator Analysis: P_3 P_1 P_2 and P_2 P_3 P_1 P_2
where P_j = (1 + s_j)/2 are symmetrizers for the j-th simple transposition.

We work in Young's seminormal form for S_4.
The irreps are indexed by partitions of 4:
  (4)      — trivial, dim 1
  (1,1,1,1) — sign, dim 1
  (3,1)    — standard, dim 3
  (2,1,1)  — standard ⊗ sign, dim 3
  (2,2)    — dim 2

In seminormal form, the generator s_i acts on standard tableaux basis vectors.
The matrix entry for s_i on tableau T is:
  s_i · T = (1/d) T + sqrt(1 - 1/d^2) T'
where d = content(i+1) - content(i) is the axial distance,
and T' is obtained by swapping i and i+1 in T (if it is standard).

Content of a box in row r, column c (0-indexed): c - r.
"""

import numpy as np
from fractions import Fraction
from itertools import permutations
import sympy as sp

# ─── Standard Young Tableaux for each partition of 4 ───────────────────────

def content(tableau, val):
    """Content (c - r, 0-indexed) of the box containing val in tableau."""
    for r, row in enumerate(tableau):
        for c, v in enumerate(row):
            if v == val:
                return c - r
    raise ValueError(f"{val} not in tableau")

def swap_values(tableau, a, b):
    """Return new tableau with a and b swapped."""
    result = []
    for row in tableau:
        new_row = []
        for v in row:
            if v == a:
                new_row.append(b)
            elif v == b:
                new_row.append(a)
            else:
                new_row.append(v)
        result.append(tuple(new_row))
    return tuple(result)

def is_standard(tableau):
    """Check if tableau is standard (rows and columns strictly increasing)."""
    for row in tableau:
        for i in range(len(row) - 1):
            if row[i] >= row[i+1]:
                return False
    n_rows = len(tableau)
    max_cols = max(len(row) for row in tableau)
    for c in range(max_cols):
        col_vals = [tableau[r][c] for r in range(n_rows) if c < len(tableau[r])]
        for i in range(len(col_vals) - 1):
            if col_vals[i] >= col_vals[i+1]:
                return False
    return True

def standard_tableaux(shape):
    """Generate all standard Young tableaux of given shape."""
    n = sum(shape)
    results = []

    def fill(partial, remaining, pos):
        """Fill tableau recursively."""
        if not remaining:
            results.append(tuple(tuple(row) for row in partial))
            return

        # Try placing next value in each valid position
        val = n - len(remaining) + 1
        for r in range(len(shape)):
            c = sum(1 for x in partial[r] if x > 0)
            if c >= shape[r]:
                continue
            # Check: this cell is at end of its row so far,
            # and either top row or cell above is filled
            if r > 0 and len(partial[r]) >= len(partial[r-1]):
                continue
            # place val here
            partial[r].append(val)
            fill(partial, remaining[1:], pos + 1)
            partial[r].pop()

    partial = [[] for _ in shape]
    fill(partial, list(range(1, n+1)), 0)
    return results

def seminormal_generator(shape, i, tableaux):
    """
    Compute matrix of s_i in Young's seminormal form.

    For each tableau T in the basis:
    - d = content(i+1) - content(i) in T
    - s_i · T = (1/d) T + sqrt(1 - 1/d^2) T_{swap}
      (where T_{swap} is T with i and i+1 swapped, if standard)
    - If T_{swap} is not standard: s_i · T = ±T (eigenvalue +1 or -1)
      specifically +1 if d > 0 (i+1 is to the right of i),
               but wait: need to be careful.
      Actually: if swapping i,i+1 doesn't give a standard tableau,
      then T is an eigenvector of s_i with eigenvalue sign(d) ...
      No: eigenvalue is +1 if i+1 is directly right (d=1) — impossible since swap would be standard

    Let me redo: the seminormal matrix entries are:
      <T | s_i | T>   = 1/d(i,T)
      <T'| s_i | T>   = sqrt(1 - 1/d^2)  if T' = swap(T,i,i+1) is standard
    where d(i,T) = content_T(i+1) - content_T(i).

    Special cases: if |d|=1 the off-diagonal term is 0 (s_i fixes T up to sign).
    Wait: sqrt(1 - 1/1) = 0, correct.

    The formula works for all d ≠ 0.
    (d=0 can't happen for standard tableaux.)
    """
    idx = {T: k for k, T in enumerate(tableaux)}
    dim = len(tableaux)
    mat = sp.zeros(dim, dim)

    for k, T in enumerate(tableaux):
        d = content(T, i+1) - content(T, i)
        # Diagonal entry
        mat[k, k] = sp.Rational(1, d)

        # Off-diagonal: try swap
        T_swap = swap_values(T, i, i+1)
        if T_swap != T and is_standard(T_swap):
            j = idx[T_swap]
            off = sp.sqrt(1 - sp.Rational(1, d*d))
            mat[k, j] = off

    return mat

# ─── Build all irreps ────────────────────────────────────────────────────────

shapes = {
    '(4)':       (4,),
    '(3,1)':     (3, 1),
    '(2,2)':     (2, 2),
    '(2,1,1)':   (2, 1, 1),
    '(1,1,1,1)': (1, 1, 1, 1),
}

def build_irrep(shape):
    tabs = standard_tableaux(shape)
    s = {}
    for i in range(1, 4):  # s_1, s_2, s_3
        s[i] = seminormal_generator(shape, i, tabs)
    return s, tabs

def proj(s_i_mat):
    """P_i = (1 + s_i) / 2"""
    dim = s_i_mat.shape[0]
    return (sp.eye(dim) + s_i_mat) / 2

def analyze_operator(name, M, s1_mat, s3_mat, label=""):
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    print(f"Matrix ({label}):")
    sp.pprint(M)

    # Rank and kernel
    rank = M.rank()
    ker = M.nullspace()
    print(f"\nRank: {rank}")
    print(f"Nullity: {len(ker)}")
    if ker:
        print("Kernel basis:")
        for v in ker:
            sp.pprint(v.T)

    # Eigenvalues
    evs = M.eigenvals()
    print(f"Eigenvalues: {evs}")

    # Check: is kernel ⊂ -1 eigenspace of s_1?
    # i.e., for each ker vector v, is s_1 v = -v?
    print("\nKernel analysis:")
    if ker:
        # eigenspace of s_1 for eigenvalue -1
        s1_minus_eig = (s1_mat + sp.eye(s1_mat.shape[0])).nullspace()  # s_1 v = -v → (s_1+I)v=0
        # Check each ker vector
        ker_in_s1_minus = True
        for v in ker:
            # Is (s_1 + I) v = 0?
            check = (s1_mat + sp.eye(s1_mat.shape[0])) * v
            if check != sp.zeros(check.shape[0], 1):
                ker_in_s1_minus = False
                break
        print(f"  kernel ⊂ (-1)-eigenspace of s_1: {ker_in_s1_minus}")
    else:
        print("  (trivial kernel)")

    # Check: is image ⊂ +1 eigenspace of s_3 AND +1 eigenspace of s_1?
    # image = column space of M
    img_basis = M.columnspace()
    print("\nImage analysis:")
    if img_basis:
        img_in_s3_plus = True
        img_in_s1_plus = True
        for v in img_basis:
            # s_3 v = +v?
            check3 = s3_mat * v - v
            if check3 != sp.zeros(check3.shape[0], 1):
                img_in_s3_plus = False
            # s_1 v = +v?
            check1 = s1_mat * v - v
            if check1 != sp.zeros(check1.shape[0], 1):
                img_in_s1_plus = False
        print(f"  image ⊂ (+1)-eigenspace of s_3: {img_in_s3_plus}")
        print(f"  image ⊂ (+1)-eigenspace of s_1: {img_in_s1_plus}")
        print(f"  image ⊂ (+1)(s_3) ∩ (+1)(s_1): {img_in_s3_plus and img_in_s1_plus}")
    else:
        print("  (zero image)")

    return rank, ker, evs

# ─── Main computation ────────────────────────────────────────────────────────

print("S_4 OPERATOR ANALYSIS")
print("Computing P_3 P_1 P_2 and P_2 P_3 P_1 P_2")
print("P_j = (1 + s_j)/2")
print()

results = {}

for shape_name, shape in shapes.items():
    print(f"\n{'#'*60}")
    print(f"# IRREP {shape_name}  (dim {sum(1 for _ in standard_tableaux(shape))})")
    print(f"{'#'*60}")

    s, tabs = build_irrep(shape)

    print(f"\nStandard tableaux ({len(tabs)}):")
    for k, T in enumerate(tabs):
        print(f"  T_{k+1} = {T}")

    # Print generator matrices
    for i in [1, 2, 3]:
        print(f"\ns_{i}:")
        sp.pprint(s[i])

    # Projectors
    P = {i: proj(s[i]) for i in [1, 2, 3]}

    # P_3 P_1 P_2
    op1 = P[3] * P[1] * P[2]

    # P_2 P_3 P_1 P_2
    op2 = P[2] * P[3] * P[1] * P[2]

    print(f"\n--- P_3 P_1 P_2 ---")
    r1, k1, e1 = analyze_operator(
        f"P_3 P_1 P_2 in irrep {shape_name}",
        op1, s[1], s[3], shape_name
    )

    print(f"\n--- P_2 P_3 P_1 P_2 ---")
    r2, k2, e2 = analyze_operator(
        f"P_2 P_3 P_1 P_2 in irrep {shape_name}",
        op2, s[1], s[3], shape_name
    )

    results[shape_name] = {
        'dim': len(tabs),
        'P3P1P2': {'rank': r1, 'ker': k1, 'evs': e1},
        'P2P3P1P2': {'rank': r2, 'ker': k2, 'evs': e2},
    }

# ─── Summary table ───────────────────────────────────────────────────────────

print("\n\n" + "="*70)
print("SUMMARY TABLE")
print("="*70)
print(f"{'Irrep':<12} {'dim':<5} {'rk(P3P1P2)':<12} {'evs(P3P1P2)':<30} {'rk(P2P3P1P2)':<15}")
print("-"*70)
for name, res in results.items():
    r1 = res['P3P1P2']['rank']
    r2 = res['P2P3P1P2']['rank']
    e1 = res['P3P1P2']['evs']
    print(f"{name:<12} {res['dim']:<5} {r1:<12} {str(e1):<30} {r2:<15}")
