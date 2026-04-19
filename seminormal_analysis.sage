"""
Seminormal form analysis of the staircase product image.

Goal: Understand which SYT basis vectors survive in im(Pi_{k-1})
and why the image is transverse to -1(s_k) eigenvectors for even k >= 4.

In Young's seminormal form:
  s_j v_T = (1/rho) v_T + sqrt(1 - 1/rho^2) v_{s_j.T}
where rho = c_T(j+1) - c_T(j) = axial distance.
If s_j.T is not standard, then s_j v_T = (1/rho) v_T with rho = +/-1.
"""

from sage.all import *
import numpy as np

def axial_distance(T, k):
    """Axial distance rho(T,k) = content(pos of k+1) - content(pos of k).
    content of cell (r,c) = c - r (0-indexed)."""
    # Find positions of k and k+1 in T
    n = T.size()
    pos_k = None
    pos_k1 = None
    for i, row in enumerate(T):
        for j, val in enumerate(row):
            if val == k:
                pos_k = (i, j)
            elif val == k + 1:
                pos_k1 = (i, j)
    ck = pos_k[1] - pos_k[0]
    ck1 = pos_k1[1] - pos_k1[0]
    return ck1 - ck

def seminormal_matrix(lam, k):
    """Build the seminormal form matrix for s_k acting on V_lambda.
    Basis = list of SYT of shape lambda, in standard order."""
    tabs = StandardTableaux(lam).list()
    d = len(tabs)
    tab_index = {T: i for i, T in enumerate(tabs)}
    
    M = np.zeros((d, d))
    for T in tabs:
        i = tab_index[T]
        rho = axial_distance(T, k)
        M[i, i] = 1.0 / rho
        
        # Try swapping k and k+1
        T_list = [list(row) for row in T]
        # Find positions
        pos_k = pos_k1 = None
        for r, row in enumerate(T_list):
            for c, val in enumerate(row):
                if val == k:
                    pos_k = (r, c)
                elif val == k + 1:
                    pos_k1 = (r, c)
        
        # Swap
        T_swap = [list(row) for row in T_list]
        T_swap[pos_k[0]][pos_k[1]] = k + 1
        T_swap[pos_k1[0]][pos_k1[1]] = k
        
        try:
            T_new = StandardTableau([tuple(row) for row in T_swap])
            # Check if it's actually standard
            if T_new in tab_index:
                j = tab_index[T_new]
                coeff = float(np.sqrt(1.0 - 1.0/rho**2))
                M[i, j] = coeff
        except:
            pass
    
    return M, tabs

def projection_matrix(lam, k):
    """P_k = (1 + s_k)/2"""
    M, tabs = seminormal_matrix(lam, k)
    d = len(tabs)
    return (np.eye(d) + M) / 2.0, tabs

def staircase_product(lam, up_to_block=None):
    """Compute staircase product Pi for S_n on V_lambda.
    Blocks: block k = (s_k, s_{k-1}, ..., s_1) for k = 1, ..., n-1.
    Product applied left-to-right in the staircase word."""
    n = sum(lam)
    if up_to_block is None:
        up_to_block = n - 1
    
    d = StandardTableaux(lam).cardinality()
    result = np.eye(d)
    tabs = None
    
    for block in range(1, up_to_block + 1):
        for k in range(block, 0, -1):
            P, tabs = projection_matrix(lam, k)
            result = P @ result
    
    return result, tabs

def minus1_eigenspace(lam, k):
    """Compute -1 eigenspace of s_k on V_lambda.
    Returns basis vectors as columns of a matrix."""
    M, tabs = seminormal_matrix(lam, k)
    d = len(tabs)
    # -1 eigenspace = null space of (s_k + I)
    A = M + np.eye(d)
    U, s, Vt = np.linalg.svd(A)
    # Null space = rows of Vt corresponding to zero singular values
    null_mask = s < 1e-10
    null_basis = Vt[null_mask].T  # columns are null vectors... wait
    # Actually null space of A is spanned by rows of Vt with zero singular values
    # But we need ker(A) = {x : Ax = 0}
    # SVD: A = U S Vt, so Ax = 0 iff S Vt x = 0 iff Vt x has zero in non-null positions
    # ker(A) = rows of Vt with singular value 0, transposed to columns
    null_basis = Vt[null_mask].T
    return null_basis, tabs

# ========== Analysis for n=5, block 4 (the first gap case) ==========
print("=" * 70)
print("SEMINORMAL FORM ANALYSIS: n=5, block k=4 (even, first gap)")
print("=" * 70)

for lam in Partitions(5):
    lam = list(lam)
    tabs = StandardTableaux(lam).list()
    d = len(tabs)
    if d <= 1:
        continue
    
    # Compute Pi_3 (staircase through block 3)
    Pi3, tabs = staircase_product(lam, up_to_block=3)
    
    # Image of Pi_3
    U, s, Vt = np.linalg.svd(Pi3)
    rank_mask = s > 1e-10
    im_Pi3 = U[:, :len(s)][:, rank_mask]  # columns = basis of image
    rank_Pi3 = im_Pi3.shape[1]
    
    # -1 eigenspace of s_4
    minus1_basis, _ = minus1_eigenspace(lam, 4)
    dim_minus1 = minus1_basis.shape[1] if len(minus1_basis.shape) > 1 else 0
    
    if dim_minus1 == 0:
        print(f"\nlambda = {lam}: dim V = {d}, rk(Pi_3) = {rank_Pi3}, dim(-1(s4)) = 0 [trivial]")
        continue
    
    # Intersection: im(Pi_3) ∩ ker(s_4 + 1)
    # A vector in both subspaces: v = im_Pi3 @ a = minus1_basis @ b
    # => [im_Pi3 | -minus1_basis] [a; b] = 0
    combined = np.hstack([im_Pi3, -minus1_basis])
    _, s_comb, Vt_comb = np.linalg.svd(combined)
    null_dim = max(0, combined.shape[1] - np.sum(s_comb > 1e-10))
    
    print(f"\nlambda = {lam}: dim V = {d}, rk(Pi_3) = {rank_Pi3}, "
          f"dim(-1(s4)) = {dim_minus1}, intersection dim = {null_dim}")
    
    # Which SYT basis vectors have nonzero coefficients in im(Pi_3)?
    print(f"  Image basis coefficients in SYT basis:")
    for col in range(rank_Pi3):
        v = im_Pi3[:, col]
        nonzero = [(i, v[i]) for i in range(d) if abs(v[i]) > 1e-10]
        tab_labels = [str(tabs[i]) for i, c in nonzero]
        print(f"    basis vec {col}: support on {len(nonzero)} SYT")
    
    # -1 eigenvectors of s_4: which SYT do they involve?
    print(f"  -1(s4) eigenvectors:")
    for col in range(dim_minus1):
        v = minus1_basis[:, col]
        nonzero = [(i, v[i]) for i in range(d) if abs(v[i]) > 1e-10]
        print(f"    -1 eigvec {col}: support on SYT indices {[i for i,c in nonzero]}")
        for i, c in nonzero:
            print(f"      T={tabs[i]}, coeff={c:.6f}, axial_dist(T,4)={axial_distance(tabs[i],4)}")

