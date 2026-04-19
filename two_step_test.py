"""
Refined gap analysis: the cascade transport works for j <= k-3 by 
commutativity. The only problematic step is j = k-2 (first cascade step).

But the image at that point is P_{k-1}(W) where W entered block k-1
already in +1(s_1).

Two-step condition: P_{k-2}P_{k-1}(+1(s_1)) ∩ -1(s_k) = {0} ?

If this holds, the gap is CLOSED because:
- After P_{k-1}: no -1(s_k) by S3 [proved]
- After P_{k-2}: the two-step condition gives no -1(s_k) [if this holds]
- After P_j for j <= k-3: commutativity preserves no -1(s_k) [proved]
"""
import numpy as np

# ========== SYT and seminormal form infrastructure ==========
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
    d = M.shape[0]
    A = M - eigenval * np.eye(d)
    U, s, Vt = np.linalg.svd(A)
    null_mask = s < tol
    if not np.any(null_mask):
        return np.zeros((d, 0))
    return Vt[null_mask].T

def image_basis(M, tol=1e-10):
    U, s, Vt = np.linalg.svd(M, full_matrices=True)
    rank = np.sum(s > tol)
    if rank == 0:
        return np.zeros((M.shape[0], 0))
    return U[:, :rank]

def subspace_intersection_dim(A, B, tol=1e-8):
    if A.shape[1] == 0 or B.shape[1] == 0:
        return 0
    combined = np.hstack([A, -B])
    _, s, _ = np.linalg.svd(combined)
    return max(0, combined.shape[1] - np.sum(s > tol))

# ============================================================
# TEST: P_{k-2}P_{k-1}(+1(s_1)) ∩ -1(s_k) = {0} ?
# For k = 4: P_2 P_3 (+1(s_1)) ∩ -1(s_4) = {0} ?
# For k = 6: P_4 P_5 (+1(s_1)) ∩ -1(s_6) = {0} ?
# ============================================================

print("=" * 70)
print("TWO-STEP CONDITION: P_{k-2}P_{k-1}(+1(s_1)) ∩ -1(s_k) = {0} ?")
print("=" * 70)

for n in range(5, 9):
    print(f"\n--- n = {n} ---")
    for lam in partitions(n):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        if d <= 1:
            continue
        
        for k in range(4, n, 2):  # even k: 4, 6, ...
            if k >= n:
                continue
            
            S1, _ = seminormal_matrix(lam, 1)
            Sk_2, _ = seminormal_matrix(lam, k-2)
            Sk_1, _ = seminormal_matrix(lam, k-1)
            Sk, _ = seminormal_matrix(lam, k)
            
            P_k2 = (np.eye(d) + Sk_2) / 2.0
            P_k1 = (np.eye(d) + Sk_1) / 2.0
            
            # +1(s_1)
            plus1_s1 = eigenspace(S1, 1.0)
            if plus1_s1.shape[1] == 0:
                continue
            
            # P_{k-1}(+1(s_1))
            img1 = P_k1 @ plus1_s1
            img1_basis = image_basis(img1)
            
            # P_{k-2}P_{k-1}(+1(s_1))
            img2 = P_k2 @ img1_basis
            img2_basis = image_basis(img2)
            
            if img2_basis.shape[1] == 0:
                continue
            
            # -1(s_k)
            minus1_sk = eigenspace(Sk, -1.0)
            if minus1_sk.shape[1] == 0:
                continue
            
            inter = subspace_intersection_dim(img2_basis, minus1_sk)
            
            status = "*** FAIL ***" if inter > 0 else "[OK]"
            if inter > 0 or n <= 6:
                print(f"  lambda={lam}, k={k}: "
                      f"P_{k-2}P_{k-1}(+1(s_1)) dim={img2_basis.shape[1]}, "
                      f"-1(s_{k}) dim={minus1_sk.shape[1]}, "
                      f"inter={inter} {status}")

print("\n" + "=" * 70)
print("If all [OK], the TWO-STEP condition holds generically on +1(s_1).")
print("Combined with commutativity transport for j <= k-3, this CLOSES THE GAP.")
print("=" * 70)

# ============================================================
# Also test the FULL cascade: P_1...P_{k-2}P_{k-1}(+1(s_1)) ∩ -1(s_k) = {0}
# This should also work and is the actual condition at block start.
# ============================================================
print("\n\nFULL CASCADE TEST: P_1...P_{k-2}P_{k-1}(+1(s_1)) ∩ -1(s_k) = {0}")
print("=" * 70)

for n in range(5, 9):
    print(f"\n--- n = {n} ---")
    for lam in partitions(n):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        if d <= 1:
            continue
        
        for k in range(4, n, 2):
            if k >= n:
                continue
            
            S1, _ = seminormal_matrix(lam, 1)
            Sk, _ = seminormal_matrix(lam, k)
            
            plus1_s1 = eigenspace(S1, 1.0)
            if plus1_s1.shape[1] == 0:
                continue
            
            # Full cascade P_1...P_{k-2}P_{k-1} applied to +1(s_1)
            W = plus1_s1.copy()
            for j in range(k-1, 0, -1):
                Sj, _ = seminormal_matrix(lam, j)
                Pj = (np.eye(d) + Sj) / 2.0
                W = Pj @ W
            W_basis = image_basis(W)
            
            if W_basis.shape[1] == 0:
                continue
            
            minus1_sk = eigenspace(Sk, -1.0)
            if minus1_sk.shape[1] == 0:
                continue
            
            inter = subspace_intersection_dim(W_basis, minus1_sk)
            status = "*** FAIL ***" if inter > 0 else "[OK]"
            if inter > 0 or n <= 6:
                print(f"  lambda={lam}, k={k}: cascade dim={W_basis.shape[1]}, "
                      f"-1(s_{k}) dim={minus1_sk.shape[1]}, inter={inter} {status}")

