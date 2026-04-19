"""
Verify Conjecture 4 for k=12 on the SKIPPED large S_13 irreps.
Use iterative image tracking: maintain an orthonormal basis for im(Pi)
and update it step by step, avoiding full d×d products.
"""
import numpy as np
import time

def partitions(n):
    if n == 0: yield (); return
    def _gen(n, mv):
        if n == 0: yield (); return
        for i in range(min(n, mv), 0, -1):
            for r in _gen(n-i, i): yield (i,) + r
    yield from _gen(n, n)

def standard_tableaux(lam):
    n = sum(lam)
    if n == 0: return [[]]
    results = []
    def fill(t, v):
        if v > n: results.append([list(r) for r in t]); return
        for i in range(len(lam)):
            if len(t[i]) < lam[i]:
                if i == 0 or len(t[i]) < len(t[i-1]):
                    t[i].append(v); fill(t, v+1); t[i].pop()
    fill([[] for _ in range(len(lam))], 1)
    return results

def find_entry(T, val):
    for r, row in enumerate(T):
        for c, v in enumerate(row):
            if v == val: return (r, c)

def axial_distance(T, k):
    r1, c1 = find_entry(T, k); r2, c2 = find_entry(T, k+1)
    return (c2 - r2) - (c1 - r1)

def tab_key(T): return tuple(tuple(r) for r in T)

def swap_entries(T, k):
    T_new = [list(r) for r in T]
    r1, c1 = find_entry(T, k); r2, c2 = find_entry(T, k+1)
    T_new[r1][c1] = k+1; T_new[r2][c2] = k
    for r, row in enumerate(T_new):
        for c, v in enumerate(row):
            if c > 0 and v <= row[c-1]: return None
            if r > 0 and c < len(T_new[r-1]) and v <= T_new[r-1][c]: return None
    return T_new

def build_sparse_seminormal(lam, k, tabs, tab_idx):
    """Return (diag, pairs) representation of s_k.
    diag[i] = 1/rho_i for fixed-point SYTs.
    pairs = list of (i, j, rho) for paired SYTs."""
    d = len(tabs)
    diag = np.zeros(d)
    pairs = []
    seen = set()
    
    for T in tabs:
        i = tab_idx[tab_key(T)]
        if i in seen:
            continue
        rho = axial_distance(T, k)
        diag[i] = 1.0 / rho
        
        T_sw = swap_entries(T, k)
        if T_sw is not None:
            key = tab_key(T_sw)
            if key in tab_idx:
                j = tab_idx[key]
                pairs.append((i, j, rho))
                seen.add(i)
                seen.add(j)
    
    return diag, pairs

def apply_projection(Q, diag, pairs, d):
    """Apply P_k = (1+s_k)/2 to the subspace spanned by columns of Q.
    Returns new Q (orthonormal basis of image).
    
    s_k acts as: s_k e_i = (1/rho) e_i + sqrt(1-1/rho^2) e_j  for pair (i,j)
                 s_k e_i = (1/rho) e_i  for fixed points
    
    P_k = (1+s_k)/2
    """
    r = Q.shape[1]  # current rank
    if r == 0:
        return Q
    
    # Compute P_k @ Q column by column
    # P_k v = (v + s_k v) / 2
    # s_k v[i] = diag[i] * v[i] + off-diag contributions
    
    result = np.zeros_like(Q)
    
    for col in range(r):
        v = Q[:, col]
        sv = np.zeros(d)
        
        # Apply s_k to v
        # Start with diagonal part
        sv = diag * v
        
        # Add off-diagonal contributions from pairs
        for (i, j, rho) in pairs:
            c = np.sqrt(1.0 - 1.0/rho**2)
            sv[i] = (1.0/rho) * v[i] + c * v[j]
            sv[j] = c * v[i] + (-1.0/rho) * v[j]
        
        result[:, col] = (v + sv) / 2.0
    
    # Re-orthogonalize and find rank
    Q_new, R = np.linalg.qr(result, mode='reduced')
    rank_new = np.sum(np.abs(np.diag(R)) > 1e-10)
    
    return Q_new[:, :rank_new]

def verify_large_irrep(lam, k):
    """Verify Conjecture 4 for a single large irrep."""
    tabs = standard_tableaux(lam)
    d = len(tabs)
    tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}
    
    # Start with full space (identity = all columns)
    Q = np.eye(d)  # This is d x d, still large...
    
    # Actually, let's track the image iteratively
    # Start: Q = I (full space)
    # After each P_j: Q = orth_basis(P_j @ Q)
    
    # The key insight: the rank drops quickly, so Q becomes thin
    # After block 1 (P_1): rank = dim(+1(s_1)) ≈ d/2
    # After block 3 (P_3, H-gen): rank drops further
    # By block k-1, rank = dim(V^H_{k-1}) which is typically small
    
    # So we need to track through the staircase, maintaining Q as thin as possible
    
    # Build all seminormal data
    sem_data = {}
    for j in range(1, k+1):
        if j < len(lam) + max(lam):  # generator j needs positions j, j+1
            sem_data[j] = build_sparse_seminormal(lam, j, tabs, tab_idx)
    
    # Apply staircase block by block
    Q = np.eye(d)
    
    for block in range(1, k):
        for j in range(block, 0, -1):
            if j not in sem_data:
                continue
            diag, pairs = sem_data[j]
            Q = apply_projection(Q, diag, pairs, d)
            if Q.shape[1] == 0:
                return True, 0, 0  # rank 0, trivially ok
    
    rk = Q.shape[1]
    
    # Now compute -1(s_k) eigenspace
    if k not in sem_data:
        return True, rk, 0
    
    diag_k, pairs_k = sem_data[k]
    
    # -1 eigenspace of s_k
    # For fixed points with rho = -1: e_i is a -1 eigenvector
    # For pairs (i,j,rho): -1 eigenvector is proportional to 
    #   -sqrt(1-1/rho^2) e_i + (1/rho + 1) e_j (... need to compute properly)
    
    # Actually, let me use full matrix for s_k eigenspace
    M_k = np.zeros((d, d))
    M_k[np.arange(d), np.arange(d)] = diag_k
    for (i, j, rho) in pairs_k:
        c = np.sqrt(1.0 - 1.0/rho**2)
        M_k[i, i] = 1.0/rho
        M_k[i, j] = c
        M_k[j, i] = c
        M_k[j, j] = -1.0/rho
    
    A = M_k + np.eye(d)
    U, s, Vt = np.linalg.svd(A)
    null_mask = s < 1e-10
    if not np.any(null_mask):
        return True, rk, 0
    m1_basis = Vt[null_mask].T
    dim_m1 = m1_basis.shape[1]
    
    # Intersection of im(Q) and -1(s_k)
    combined = np.hstack([Q, -m1_basis])
    _, sc, _ = np.linalg.svd(combined)
    inter = max(0, combined.shape[1] - np.sum(sc > 1e-8))
    
    return inter == 0, rk, dim_m1

# Test on skipped partitions
k = 12
n = k + 1
skipped = []
for lam in partitions(n):
    d = len(standard_tableaux(lam))
    if d > 5000 and d <= 1:
        continue
    if d > 5000:
        skipped.append((lam, d))

print(f"Skipped partitions for k={k} (S_{n}): {len(skipped)}")
print(f"Dimensions: {[d for _, d in skipped]}")
print(f"Will attempt iterative verification...")
print()

# Try smallest skipped first
skipped.sort(key=lambda x: x[1])

for lam, d in skipped[:5]:  # try first 5 smallest
    t0 = time.time()
    print(f"  lambda={lam}, dim={d}...", end=" ", flush=True)
    try:
        ok, rk, dim_m1 = verify_large_irrep(lam, k)
        t1 = time.time()
        status = "OK" if ok else "FAIL"
        print(f"rk={rk}, dim(-1)={dim_m1}, [{status}] ({t1-t0:.1f}s)")
    except Exception as e:
        t1 = time.time()
        print(f"ERROR: {e} ({t1-t0:.1f}s)")

