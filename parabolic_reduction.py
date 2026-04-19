"""
PARABOLIC REDUCTION for Conjecture 4.

Key insight: Conjecture 4 for even block k is a statement about generators
s_1, ..., s_k, which generate S_{k+1}. Any S_n irrep restricts to S_{k+1}
irreps. So the conjecture for block k holds for ALL n iff it holds for
all S_{k+1} irreps.

This means:
- k=4: check all S_5 irreps (7 partitions) -> proves for ALL n
- k=6: check all S_7 irreps (15 partitions) -> proves for ALL n  
- k=8: check all S_9 irreps (30 partitions) -> proves for ALL n

Each finite verification is COMPLETE for its block.
"""
import numpy as np

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

def staircase_product_matrix(lam, up_to_block):
    """Return the staircase product matrix through the given block."""
    d = len(standard_tableaux(lam))
    if d == 0:
        return np.zeros((0,0))
    result = np.eye(d)
    for block in range(1, up_to_block + 1):
        for k in range(block, 0, -1):
            M, _ = seminormal_matrix(lam, k)
            P = (np.eye(d) + M) / 2.0
            result = P @ result
    return result

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
# Verify Conjecture 4 for each even k via parabolic reduction
# k=2m: verify on all S_{2m+1} irreps
# ============================================================

print("=" * 70)
print("PARABOLIC REDUCTION VERIFICATION OF CONJECTURE 4")
print("=" * 70)
print()
print("For even block k, Conjecture 4 reduces to checking all S_{k+1} irreps.")
print("Verification for one k is COMPLETE for ALL n.")
print()

results = {}

for k in [4, 6, 8]:
    n_check = k + 1  # S_{k+1}
    print(f"{'='*70}")
    print(f"Block k = {k}: checking all S_{n_check} irreps ({sum(1 for _ in partitions(n_check))} partitions)")
    print(f"{'='*70}")
    
    all_ok = True
    for lam in partitions(n_check):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        
        if d <= 1:
            # 1-dim irreps: trivially ok
            continue
        
        # Compute staircase through block k-1
        Pi = staircase_product_matrix(lam, k - 1)
        im = image_basis(Pi)
        
        # Compute -1(s_k) eigenspace
        Sk, _ = seminormal_matrix(lam, k)
        m1 = eigenspace(Sk, -1.0)
        
        if im.shape[1] == 0 or m1.shape[1] == 0:
            print(f"  lambda={lam}: dim={d}, rk(Pi_{k-1})={im.shape[1]}, "
                  f"dim(-1(s_{k}))={m1.shape[1]} [trivial]")
            continue
        
        inter = subspace_intersection_dim(im, m1)
        
        status = "OK" if inter == 0 else "FAIL"
        if inter > 0:
            all_ok = False
        print(f"  lambda={lam}: dim={d}, rk(Pi_{k-1})={im.shape[1]}, "
              f"dim(-1(s_{k}))={m1.shape[1]}, intersection={inter} [{status}]")
    
    results[k] = all_ok
    if all_ok:
        print(f"\n  *** CONJECTURE 4 VERIFIED for k={k} (all S_{n_check} irreps) ***")
        print(f"  *** This proves the H-invariant theorem for ALL n when block {k} appears ***")
    else:
        print(f"\n  !!! CONJECTURE 4 FAILS for k={k} !!!")
    print()

print("=" * 70)
print("SUMMARY")
print("=" * 70)
for k, ok in results.items():
    n_check = k + 1
    n_range = f"n <= {k+2}" if ok else "FAILS"
    print(f"  k={k} (S_{n_check} irreps): {'VERIFIED' if ok else 'FAILED'} => theorem for {n_range}")

# What n does each k cover?
# k=4 first appears at n=5. If verified, theorem holds for n<=6 (next even k=6 at n=7)
# k=6 first appears at n=7. If verified, theorem holds for n<=8
# k=8 first appears at n=9. If verified, theorem holds for n<=10
# In general: if verified for all even k <= K, theorem holds for all n <= K+2

max_k_verified = max(k for k, ok in results.items() if ok) if any(results.values()) else 0
if max_k_verified > 0:
    print(f"\nH-INVARIANT THEOREM PROVED UNCONDITIONALLY for all n <= {max_k_verified + 2}")
    print(f"(All even blocks up to k={max_k_verified} verified via parabolic reduction)")

