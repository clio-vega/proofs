"""Compute dim B^{(3,2,1)} at indices 4 (sharp), 3 (one below), and check E_1^- containment."""
import sys
sys.path.insert(0, '/home/clio/projects/scratch/2026-05-13-multi-corner')
import numpy as np
from compute_dim_B import (
    build_T_i_p1_modp, build_T_i_modp, syt_of_shape,
    matrix_rank_modp, matmul_modp_safe, P
)

def compute_B_image(lam, n, j0_target):
    syts = syt_of_shape(lam)
    ndim = len(syts)
    Ti_p1 = {}
    for i in range(1, n):
        Ti_p1[i] = build_T_i_p1_modp(syts, i)
    M = np.eye(ndim, dtype=np.int64)
    for k in range(n, j0_target - 1, -1):
        Rk = np.eye(ndim, dtype=np.int64)
        for i in range(1, k):
            Rk = matmul_modp_safe(Ti_p1[i], Rk)
        M = matmul_modp_safe(Rk, M)
    return M  # image-rank columns

lam = (3, 2, 1)
n = 6
syts = syt_of_shape(lam)
ndim = len(syts)
print(f"V_{lam} at n={n}: dim {ndim}")

# Compute B at various indices
for j in [4, 3, 2]:
    M = compute_B_image(lam, n, j)
    rk = matrix_rank_modp(M)
    print(f"  dim B_{j} = {rk}")

# Check if B_4 V is in E_1^-: i.e., (T_1 + 1) B_4 V = 0?
M4 = compute_B_image(lam, n, 4)
T1_p1 = build_T_i_p1_modp(syts, 1)
M_after_T1p1 = matmul_modp_safe(T1_p1, M4)
rk_after = matrix_rank_modp(M_after_T1p1)
print(f"  rank of (T_1+1) B_4 V = {rk_after}")
print(f"  (T_1+1) B_4 V == 0  iff  rank == 0:  {'YES, B_4 V is in E_1^-' if rk_after == 0 else 'NO, B_4 V not in E_1^-'}")
