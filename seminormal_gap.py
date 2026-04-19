"""
Seminormal form analysis of the staircase product image.
Goal: understand why im(Pi_{k-1}) ∩ -1(s_k) = {0} for even k >= 4.

Pure Python/NumPy implementation (no SageMath).
"""
import numpy as np
from itertools import permutations
from functools import lru_cache

# ========== Standard Young Tableaux ==========

def partitions(n):
    """Generate all partitions of n."""
    if n == 0:
        yield ()
        return
    if n < 0:
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
    """Generate all SYT of shape lam."""
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
                # Check column condition
                if i == 0 or len(tableau[i]) < len(tableau[i-1]):
                    tableau[i].append(val)
                    fill(tableau, val + 1)
                    tableau[i].pop()
    
    fill([[] for _ in range(len(lam))], 1)
    return results

def find_entry(T, val):
    """Find (row, col) of val in tableau T (0-indexed)."""
    for r, row in enumerate(T):
        for c, v in enumerate(row):
            if v == val:
                return (r, c)
    return None

def content(r, c):
    """Content of cell (r,c) = c - r."""
    return c - r

def axial_distance(T, k):
    """rho(T,k) = content(pos of k+1) - content(pos of k)."""
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k+1)
    return content(r2, c2) - content(r1, c1)

def swap_entries(T, k):
    """Swap k and k+1 in T. Return new tableau or None if not standard."""
    T_new = [list(row) for row in T]
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k+1)
    T_new[r1][c1] = k + 1
    T_new[r2][c2] = k
    
    # Check if standard
    for r, row in enumerate(T_new):
        for c, v in enumerate(row):
            if c > 0 and v <= row[c-1]:
                return None
            if r > 0 and c < len(T_new[r-1]) and v <= T_new[r-1][c]:
                return None
    return T_new

def tab_key(T):
    """Hashable key for a tableau."""
    return tuple(tuple(row) for row in T)

# ========== Seminormal representation matrices ==========

def seminormal_matrix(lam, k):
    """s_k matrix in Young's seminormal form on V_lambda."""
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

def projection_matrix(lam, k):
    """P_k = (1 + s_k)/2."""
    M, tabs = seminormal_matrix(lam, k)
    return (np.eye(len(tabs)) + M) / 2.0, tabs

def staircase_image(lam, up_to_block):
    """Compute image of staircase product Pi through block up_to_block.
    Returns: image basis (columns), tabs list."""
    d = len(standard_tableaux(lam))
    if d == 0:
        return np.zeros((0, 0)), []
    
    result = np.eye(d)
    tabs = None
    
    for block in range(1, up_to_block + 1):
        for k in range(block, 0, -1):
            P, tabs = projection_matrix(lam, k)
            result = P @ result
    
    # Extract image basis via SVD
    U, s, Vt = np.linalg.svd(result)
    rank = np.sum(s > 1e-10)
    im_basis = U[:, :rank]
    
    return im_basis, tabs, result

def minus1_eigenspace(lam, k):
    """Compute -1 eigenspace of s_k. Returns column basis."""
    M, tabs = seminormal_matrix(lam, k)
    d = len(tabs)
    A = M + np.eye(d)
    U, s, Vt = np.linalg.svd(A)
    null_mask = s < 1e-10
    if not np.any(null_mask):
        return np.zeros((d, 0)), tabs
    return Vt[null_mask].T, tabs

def subspace_intersection_dim(A, B):
    """Dimension of intersection of column spaces of A and B."""
    if A.shape[1] == 0 or B.shape[1] == 0:
        return 0
    combined = np.hstack([A, -B])
    _, s, _ = np.linalg.svd(combined)
    return max(0, combined.shape[1] - np.sum(s > 1e-8))

# ========== Main analysis ==========

print("=" * 70)
print("SEMINORMAL FORM ANALYSIS: Conjecture 4 (even block-start injectivity)")
print("=" * 70)

for n in [5, 6, 7]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")
    
    for lam in partitions(n):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        if d <= 1:
            continue
        
        # Check each even block >= 4
        for block_k in range(4, n, 2):  # even blocks: 4, 6, ...
            if block_k >= n:
                continue
            
            im_basis, tabs, Pi_mat = staircase_image(lam, block_k - 1)
            rank = im_basis.shape[1]
            
            m1_basis, _ = minus1_eigenspace(lam, block_k)
            dim_m1 = m1_basis.shape[1]
            
            if dim_m1 == 0:
                continue
            
            inter_dim = subspace_intersection_dim(im_basis, m1_basis)
            
            print(f"\nlambda={lam}, block k={block_k}: "
                  f"dim V={d}, rk(Pi_{block_k-1})={rank}, dim(-1(s{block_k}))={dim_m1}, "
                  f"intersection={inter_dim} {'*** GAP CASE ***' if inter_dim > 0 else '[OK]'}")
            
            if inter_dim == 0 and dim_m1 > 0:
                # This is the interesting case: gap holds
                # Analyze WHY: decompose image basis in SYT basis
                print(f"  Image vectors in SYT basis:")
                for col_idx in range(min(rank, 4)):
                    v = im_basis[:, col_idx]
                    support = [(i, v[i]) for i in range(d) if abs(v[i]) > 1e-10]
                    print(f"    im_vec {col_idx}: {len(support)} nonzero entries")
                    for idx, coeff in support:
                        T = tabs[idx]
                        rho_k = axial_distance(T, block_k) if block_k < sum(lam) else 'N/A'
                        print(f"      SYT {tab_key(T)}: coeff={coeff:.6f}, "
                              f"rho(T,{block_k})={rho_k}")
                
                print(f"  -1(s{block_k}) eigenvectors:")
                for col_idx in range(dim_m1):
                    v = m1_basis[:, col_idx]
                    support = [(i, v[i]) for i in range(d) if abs(v[i]) > 1e-10]
                    print(f"    -1 eigvec {col_idx}: {len(support)} nonzero entries")
                    for idx, coeff in support:
                        T = tabs[idx]
                        rho_k = axial_distance(T, block_k) if block_k < sum(lam) else 'N/A'
                        print(f"      SYT {tab_key(T)}: coeff={coeff:.6f}, "
                              f"rho(T,{block_k})={rho_k}")

