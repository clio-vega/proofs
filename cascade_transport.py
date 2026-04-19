"""
Test the CASCADE TRANSPORT hypothesis:

After P_{k-1}, image is in +1(s_{k-1}), so no -1(s_k) vectors (S3 lemma).
The cascade P_{k-2}, ..., P_1 follows. s_k commutes with P_j for j <= k-3.

Key question: does P_{k-2}(+1(s_{k-1})) ∩ -1(s_k) = {0} ?

This is a question about THREE consecutive generators s_{k-2}, s_{k-1}, s_k,
which generate a copy of S_4. If this holds for all S_4 irreps, the rest
of the cascade follows by commutativity.
"""
import numpy as np

# Reuse the SYT infrastructure
def partitions(n):
    if n == 0:
        yield ()
        return
    def _gen(n, max_val):
        if n == 0:
            yield ()
            return
        for i in range(min(n, max_val), 0, -1):
            for rest in _gen(n - i, i):
                yield (i,) + rest
    yield from _gen(n, n)

def standard_tableaux(lam):
    n = sum(lam)
    if n == 0:
        return [[]]
    results = []
    def fill(tableau, val):
        if val > n:
            results.append([list(row) for row in tableau])
            return
        for i in range(len(lam)):
            if len(tableau[i]) < lam[i]:
                if i == 0 or len(tableau[i]) < len(tableau[i-1]):
                    tableau[i].append(val)
                    fill(tableau, val + 1)
                    tableau[i].pop()
    fill([[] for _ in range(len(lam))], 1)
    return results

def find_entry(T, val):
    for r, row in enumerate(T):
        for c, v in enumerate(row):
            if v == val:
                return (r, c)
    return None

def axial_distance(T, k):
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k+1)
    return (c2 - r2) - (c1 - r1)

def tab_key(T):
    return tuple(tuple(row) for row in T)

def swap_entries(T, k):
    T_new = [list(row) for row in T]
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k+1)
    T_new[r1][c1] = k + 1
    T_new[r2][c2] = k
    for r, row in enumerate(T_new):
        for c, v in enumerate(row):
            if c > 0 and v <= row[c-1]:
                return None
            if r > 0 and c < len(T_new[r-1]) and v <= T_new[r-1][c]:
                return None
    return T_new

def seminormal_matrix(lam, k):
    tabs = standard_tableaux(lam)
    d = len(tabs)
    tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}
    M = np.zeros((d, d))
    for T in tabs:
        i = tab_idx[tab_key(T)]
        rho = axial_distance(T, k)
        M[i, i] = 1.0 / rho
        T_sw = swap_entries(T, k)
        if T_sw is not None:
            key = tab_key(T_sw)
            if key in tab_idx:
                j = tab_idx[key]
                M[i, j] = np.sqrt(1.0 - 1.0/rho**2)
    return M, tabs

def eigenspace(M, eigenval, tol=1e-10):
    """Column basis of eigenspace of M for given eigenvalue."""
    d = M.shape[0]
    A = M - eigenval * np.eye(d)
    U, s, Vt = np.linalg.svd(A)
    null_mask = s < tol
    if not np.any(null_mask):
        return np.zeros((d, 0))
    return Vt[null_mask].T

def subspace_intersection_dim(A, B, tol=1e-8):
    if A.shape[1] == 0 or B.shape[1] == 0:
        return 0
    combined = np.hstack([A, -B])
    _, s, _ = np.linalg.svd(combined)
    return max(0, combined.shape[1] - np.sum(s > tol))

# ============================================================
# TEST: P_{k-2}(+1(s_{k-1})) ∩ -1(s_k) = {0} ?
# Using S_n representations with consecutive generators
# ============================================================

print("=" * 70)
print("TEST: P_{j}(+1(s_{j+1})) ∩ -1(s_{j+2}) = {0} ?")
print("This is the critical single-step question.")
print("=" * 70)

for n in range(4, 9):
    print(f"\n--- n = {n} ---")
    for lam in partitions(n):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        if d <= 1:
            continue
        
        # Test for all possible triples (j, j+1, j+2) where j+2 <= n-1
        for j in range(1, n-2):
            # s_j, s_{j+1}, s_{j+2} are three consecutive generators
            S_j, _ = seminormal_matrix(lam, j)
            S_j1, _ = seminormal_matrix(lam, j+1)
            S_j2, _ = seminormal_matrix(lam, j+2)
            
            P_j = (np.eye(d) + S_j) / 2.0  # projection onto +1(s_j)
            
            # +1(s_{j+1}) eigenspace
            plus1_j1 = eigenspace(S_j1, 1.0)
            
            # -1(s_{j+2}) eigenspace
            minus1_j2 = eigenspace(S_j2, -1.0)
            
            if plus1_j1.shape[1] == 0 or minus1_j2.shape[1] == 0:
                continue
            
            # Compute P_j(+1(s_{j+1})) = image of P_j restricted to +1(s_{j+1})
            img = P_j @ plus1_j1
            # Get actual image basis via SVD
            U, s, Vt = np.linalg.svd(img, full_matrices=False)
            rank = np.sum(s > 1e-10)
            if rank == 0:
                continue
            img_basis = U[:, :rank]
            
            # Intersection with -1(s_{j+2})
            inter = subspace_intersection_dim(img_basis, minus1_j2)
            
            if inter > 0:
                print(f"  COUNTEREXAMPLE! lambda={lam}, (s_{j},s_{j+1},s_{j+2}): "
                      f"P_{j}(+1(s_{j+1})) ∩ -1(s_{j+2}) has dim {inter}")
            # Only print non-trivial cases
            elif d > 2:  # skip trivial 2-dim cases
                pass  # keep output clean

print("\n" + "=" * 70)
print("DONE. If no counterexamples, the single-step condition holds.")
print("=" * 70)

# ============================================================
# If the above holds, the full argument is:
# 1. After P_{k-1}: im in +1(s_{k-1}), no -1(s_k) by S3
# 2. P_{k-2} applied: P_{k-2}(+1(s_{k-1})) ∩ -1(s_k) = {0} [THE KEY STEP]
# 3. For j <= k-3: im in +1(s_{j+1}), s_k commutes with all of these,
#    so the "no -1(s_k)" property is preserved algebraically.
# ============================================================

print("\n\nNow test the FULL cascade transport argument:")
print("After P_{k-1}: in +1(s_{k-1}), no -1(s_k) [S3].")
print("After P_{k-2}: in +1(s_{k-2}), no -1(s_k) [KEY STEP above].")
print("After P_j (j<=k-3): s_k commutes with P_j, preserves +1(s_{j+1}).")
print("Need: +1(s_{j+1}) ∩ +1(s_k) ∩ -1(s_j) = {0} for j <= k-3.")
print()

# For j <= k-3: |k - j| >= 3, |k - (j+1)| >= 2
# s_k commutes with both s_j and s_{j+1}
# v ∈ +1(s_{j+1}), v = v+ + v- (s_k eigendecomp)
# s_k commutes with s_{j+1}, so v+ ∈ +1(s_{j+1})
# P_j v ∈ -1(s_k) iff v+ ∈ -1(s_j) (since s_k commutes with P_j)
# So need: +1(s_{j+1}) ∩ +1(s_k) ∩ -1(s_j) = {0}
# By S3 for (s_j, s_{j+1}): +1(s_{j+1}) ∩ -1(s_j) = {0}
# So this is AUTOMATIC.

print("For j <= k-3:")
print("  s_k commutes with s_j and s_{j+1}")
print("  v+ (s_k = +1 component of v) is in +1(s_{j+1}) [since s_k commutes with s_{j+1}]")
print("  P_j v ∈ -1(s_k) iff v+ ∈ -1(s_j)")
print("  But v+ ∈ +1(s_{j+1}) ∩ -1(s_j) = {0} by S3 lemma")
print("  So the cascade preserves 'no -1(s_k)' for all j <= k-3. QED for these steps.")
print()
print("The ONLY non-trivial step is j = k-2, which is the KEY STEP tested above.")

