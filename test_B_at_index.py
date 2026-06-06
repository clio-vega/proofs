"""Compute dim B at non-sharp index for hook (4,1,1)."""
import sys
sys.path.insert(0, '/home/clio/projects/scratch/2026-05-13-multi-corner')
import numpy as np
from compute_dim_B import (
    build_T_i_p1_modp, syt_of_shape, matrix_rank_modp, matmul_modp_safe, P
)

def compute_dim_B_at_index(lam, n, j0_target, verbose=False):
    """Compute dim of R'_{j0_target} ... R'_n V_lambda."""
    syts = syt_of_shape(lam)
    ndim = len(syts)
    Ti_p1 = {}
    for i in range(1, n):
        Ti_p1[i] = build_T_i_p1_modp(syts, i)
    if j0_target > n:
        return ndim, ndim
    M = np.eye(ndim, dtype=np.int64)
    for k in range(n, j0_target - 1, -1):
        Rk = np.eye(ndim, dtype=np.int64)
        for i in range(1, k):
            Rk = matmul_modp_safe(Ti_p1[i], Rk)
        M = matmul_modp_safe(Rk, M)
    return matrix_rank_modp(M), ndim

# Hook (4, 1, 1) at n = 6, sharp = j0 = 4. Compute at index 3 and 4.
print("(4, 1, 1) at n=6:")
for idx in [4, 3, 2]:
    dim, d = compute_dim_B_at_index((4,1,1), 6, idx)
    print(f"  B_{idx} dim = {dim}")

# Hook (5, 1) at n = 6, sharp = j0 = 3.
print("\n(5, 1) at n=6:")
for idx in [3, 2]:
    dim, d = compute_dim_B_at_index((5,1), 6, idx)
    print(f"  B_{idx} dim = {dim}")

# (3, 2) at n = 5, sharp = 3.
print("\n(3, 2) at n=5:")
for idx in [3, 2]:
    dim, d = compute_dim_B_at_index((3,2), 5, idx)
    print(f"  B_{idx} dim = {dim}")

# (4, 2) at n = 6, sharp = 3.
print("\n(4, 2) at n=6:")
for idx in [3, 2]:
    dim, d = compute_dim_B_at_index((4,2), 6, idx)
    print(f"  B_{idx} dim = {dim}")

# (3, 1, 1) at n=5, sharp=4.
print("\n(3, 1, 1) at n=5:")
for idx in [4, 3, 2]:
    dim, d = compute_dim_B_at_index((3,1,1), 5, idx)
    print(f"  B_{idx} dim = {dim}")
