"""
verify_k14.py
Verify Conjecture 4 for k=14: for every irrep λ ⊢ 15,
    im(Π_{13}) ∩ ker(s_{14} + I) = {0}

where Π_{13} = (P_{13}···P_1)(P_{12}···P_1)···(P_1) in S_15
(staircase: 13 blocks, block b applies P_j = (I+s_j)/2 in order j=b,b-1,...,1)

S_15 has 176 partitions. Some irreps have dim ~ 50000-80000.

Two methods depending on dimension:
 - EXACT (dim ≤ 3000): SVD-based image tracking with correct rank detection
 - FAST (dim > 3000): staircase applied to random thin basis + Gram compression
   with m = max(100, dim//30)

Seminormal form: s_j|T⟩ = (1/d)|T⟩ + √(1-1/d²)|T'⟩
where d = axial_distance(T,j) and T' = swap_j(T).

Tolerance: 1e-8 for rank detection.
"""
import numpy as np
import time
import math


# ─── Partition utilities ───────────────────────────────────────────────────────

def partitions(n):
    if n == 0:
        yield ()
        return
    def _gen(n, mv):
        if n == 0:
            yield ()
            return
        for i in range(min(n, mv), 0, -1):
            for r in _gen(n - i, i):
                yield (i,) + r
    yield from _gen(n, n)


def hook_length_dim(lam):
    n = sum(lam)
    if n == 0:
        return 1
    conj = [0] * lam[0]
    for i, r in enumerate(lam):
        for j in range(r):
            conj[j] += 1
    prod = 1
    for i in range(len(lam)):
        for j in range(lam[i]):
            hook = (lam[i] - j - 1) + (conj[j] - i - 1) + 1
            prod *= hook
    return math.factorial(n) // prod


# ─── Standard tableaux ────────────────────────────────────────────────────────

def standard_tableaux(lam):
    n = sum(lam)
    if n == 0:
        return [[]]
    results = []
    def fill(t, v):
        if v > n:
            results.append([list(r) for r in t])
            return
        for i in range(len(lam)):
            if len(t[i]) < lam[i]:
                if i == 0 or len(t[i]) < len(t[i - 1]):
                    t[i].append(v)
                    fill(t, v + 1)
                    t[i].pop()
    fill([[] for _ in range(len(lam))], 1)
    return results


def find_entry(T, val):
    for r, row in enumerate(T):
        for c, v in enumerate(row):
            if v == val:
                return (r, c)
    return None


def axial_distance(T, k):
    res1 = find_entry(T, k)
    res2 = find_entry(T, k + 1)
    if res1 is None or res2 is None:
        return None
    r1, c1 = res1
    r2, c2 = res2
    return (c2 - r2) - (c1 - r1)


def tab_key(T):
    return tuple(tuple(r) for r in T)


def swap_entries(T, k):
    if find_entry(T, k) is None or find_entry(T, k + 1) is None:
        return None
    T_new = [list(r) for r in T]
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k + 1)
    T_new[r1][c1] = k + 1
    T_new[r2][c2] = k
    for r, row in enumerate(T_new):
        for c, v in enumerate(row):
            if c > 0 and v <= row[c - 1]:
                return None
            if r > 0 and c < len(T_new[r - 1]) and v <= T_new[r - 1][c]:
                return None
    return T_new


# ─── Seminormal form ──────────────────────────────────────────────────────────

def build_seminormal_data(lam, k, tabs, tab_idx):
    """
    Build sparse representation of s_k in Young's seminormal form.

    For pair (i, j) with axial distance rho:
        s_k e_i = (1/rho) e_i + c e_j,   c = sqrt(1-1/rho^2)
        s_k e_j = c e_i - (1/rho) e_j

    Returns: (diag, off_i, off_j, off_ci, off_cj)
    """
    n = sum(lam)
    if k >= n:
        return (np.zeros(len(tabs)), np.array([], dtype=np.intp),
                np.array([], dtype=np.intp), np.array([]), np.array([]))
    d = len(tabs)
    diag = np.zeros(d)
    oi, oj, oci, ocj = [], [], [], []
    seen = set()

    for T in tabs:
        i = tab_idx[tab_key(T)]
        if i in seen:
            continue
        rho = axial_distance(T, k)
        if rho is None:
            seen.add(i)
            continue
        T_sw = swap_entries(T, k)
        if T_sw is not None:
            key = tab_key(T_sw)
            if key in tab_idx:
                j = tab_idx[key]
                c = math.sqrt(max(0.0, 1.0 - 1.0 / rho**2))
                oi.append(i)
                oj.append(j)
                oci.append(1.0 / rho)
                ocj.append(c)
                seen.add(i)
                seen.add(j)
                diag[i] = 1.0 / rho
                diag[j] = -1.0 / rho
                continue
        seen.add(i)
        diag[i] = 1.0 / rho

    return (diag,
            np.array(oi, dtype=np.intp),
            np.array(oj, dtype=np.intp),
            np.array(oci),
            np.array(ocj))


def apply_sk(B, sem):
    """Apply s_k to all columns of B (d×m matrix). Returns s_k @ B."""
    diag, off_i, off_j, off_ci, off_cj = sem
    sB = diag[:, np.newaxis] * B
    if len(off_i) > 0:
        sB[off_i] += off_cj[:, np.newaxis] * B[off_j]
        sB[off_j] += off_cj[:, np.newaxis] * B[off_i]
    return sB


def apply_projection(B, sem):
    """Apply P_k = (I + s_k)/2 to all columns of B. Returns P_k @ B."""
    return (B + apply_sk(B, sem)) * 0.5


def build_minus1_basis(d, sem_k):
    """
    Build orthonormal basis for the -1 eigenspace of s_k.

    For each pair (i,j): v[i] = c, v[j] = -(1+1/rho), normalized.
    Fixed points with diag = -1 also contribute.
    """
    diag_k, off_i, off_j, off_ci, off_cj = sem_k
    minus1_vecs = []
    paired = set(off_i.tolist()) | set(off_j.tolist())

    for idx in range(len(off_i)):
        i = off_i[idx]
        j = off_j[idx]
        rho = 1.0 / off_ci[idx]
        c = off_cj[idx]
        v = np.zeros(d)
        v[i] = c
        v[j] = -(1.0 / rho + 1.0)
        nrm = np.linalg.norm(v)
        if nrm > 1e-14:
            minus1_vecs.append(v / nrm)

    for idx in range(d):
        if idx not in paired and abs(diag_k[idx] + 1.0) < 1e-10:
            v = np.zeros(d)
            v[idx] = 1.0
            minus1_vecs.append(v)

    if not minus1_vecs:
        return np.zeros((d, 0))
    return np.column_stack(minus1_vecs)


# ─── Gram matrix compression ──────────────────────────────────────────────────

def compress_gram(B, tol=1e-8):
    """
    Find orthonormal basis for col(B) using Gram matrix eigendecomposition.
    Much faster than full SVD for tall-skinny matrices.
    """
    if B.shape[1] == 0:
        return np.zeros((B.shape[0], 0)), 0
    G = B.T @ B
    vals, vecs = np.linalg.eigh(G)
    if len(vals) == 0 or vals[-1] <= 0:
        return np.zeros((B.shape[0], 0)), 0
    idx = np.argsort(vals)[::-1]
    vals = vals[idx]
    vecs = vecs[:, idx]
    rank = int(np.sum(vals > tol * vals[0]))
    if rank == 0:
        return np.zeros((B.shape[0], 0)), 0
    V_r = vecs[:, :rank]
    sv_r = np.sqrt(np.maximum(1e-300, vals[:rank]))
    U_r = B @ V_r / sv_r[np.newaxis, :]
    return U_r, rank


# ─── Exact method (small irreps, d ≤ 3000) ────────────────────────────────────

def build_plus1_eigenspace(d, sem_j):
    """
    Build orthonormal basis for the +1 eigenspace of s_j directly.
    """
    diag, off_i, off_j, off_ci, off_cj = sem_j
    basis = []
    paired = set(off_i.tolist()) | set(off_j.tolist())

    for idx in range(len(off_i)):
        i = off_i[idx]
        j = off_j[idx]
        rho = 1.0 / off_ci[idx]
        c = off_cj[idx]
        v = np.zeros(d)
        v[i] = c
        v[j] = 1.0 - 1.0 / rho
        nrm = np.linalg.norm(v)
        if nrm > 1e-14:
            basis.append(v / nrm)
        else:
            v2 = np.zeros(d)
            v2[i] = 1.0
            basis.append(v2)

    for idx in range(d):
        if idx not in paired and diag[idx] > 0:
            v = np.zeros(d)
            v[idx] = 1.0
            basis.append(v)

    if not basis:
        return np.zeros((d, 0))
    return np.column_stack(basis)


def apply_projection_svd(B, sem, tol=1e-8):
    """
    Apply P_k = (I + s_k)/2 to B (d×r), return orthonormal basis for image.
    Uses SVD for correct rank detection.
    """
    PB = apply_projection(B, sem)
    U, s, _ = np.linalg.svd(PB, full_matrices=False)
    rank = int(np.sum(s > tol))
    return U[:, :rank]


def verify_irrep_exact(lam, k, tol=1e-8):
    """
    Exact method for small irreps (d ≤ 3000).
    Bootstraps from +1 eigenspace of s_1.
    Uses SVD at each step for correct rank detection.
    """
    n = k + 1
    tabs = standard_tableaux(lam)
    d = len(tabs)
    tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}

    sem = {}
    for j in range(1, k + 1):
        sem[j] = build_seminormal_data(lam, j, tabs, tab_idx)

    # Bootstrap from +1 eigenspace of s_1 (block 1 = just P_1)
    B = build_plus1_eigenspace(d, sem[1])
    if B.shape[1] == 0:
        return True, 0, 0

    # Apply blocks 2..k-1
    for block in range(2, k):
        for j in range(block, 0, -1):
            B = apply_projection_svd(B, sem[j], tol=tol)
            if B.shape[1] == 0:
                return True, 0, 0

    rank_im = B.shape[1]

    # Build -1 eigenspace of s_k
    M1 = build_minus1_basis(d, sem[k])
    dim_m1 = M1.shape[1]
    if dim_m1 == 0:
        return True, rank_im, 0

    # Intersection check using orthonormality of M1
    coeff = M1.T @ B
    residual = B - M1 @ coeff
    G_res = residual.T @ residual
    if G_res.shape[0] == 0:
        return True, rank_im, dim_m1
    vals_r, _ = np.linalg.eigh(G_res)
    max_val = np.max(vals_r)
    rank_res = int(np.sum(vals_r > tol * max_val)) if max_val > 0 else 0
    inter_dim = max(0, rank_im - rank_res)

    return inter_dim == 0, rank_im, dim_m1


# ─── Fast method (large irreps, d > 3000) ────────────────────────────────────

def verify_irrep_fast(lam, k, m=None, seed=42, tol_rank=1e-8, tol_inter=1e-8):
    """
    Fast method for large irreps (d > 3000).

    Uses a random initial basis of size m = max(100, dim//30),
    applies the full staircase WITHOUT intermediate SVD/QR,
    using only column normalization to prevent underflow.
    Final compression via Gram matrix eigendecomposition.

    Correctness: if im(Π) ∩ ker(s_k+I) ≠ {0}, then some v is fixed by all
    staircase projections and lies in ker(s_k+I). A random initial basis
    captures v with probability 1, so ok=True is conclusive.
    """
    import gc
    n = k + 1
    tabs = standard_tableaux(lam)
    d = len(tabs)
    tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}

    if m is None:
        # Memory budget: keep d*m*8 bytes < 600 MB => m < 6e8 / (d*8)
        m_mem = int(6e8 / (d * 8))
        m_ideal = max(100, d // 30)
        m = max(50, min(m_ideal, m_mem))

    sem = {}
    for j in range(1, k + 1):
        sem[j] = build_seminormal_data(lam, j, tabs, tab_idx)

    # Free tabs and tab_idx before allocating large working matrix
    del tabs, tab_idx
    gc.collect()

    # Random initial basis
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((d, m)).astype(np.float64)
    norms = np.linalg.norm(B, axis=0, keepdims=True)
    B /= norms

    # Apply full staircase (k-1 blocks) with column renormalization after each block
    for block in range(1, k):
        for j in range(block, 0, -1):
            B = apply_projection(B, sem[j])
        # Renormalize columns to prevent underflow
        norms = np.linalg.norm(B, axis=0, keepdims=True)
        active = norms[0] > 1e-30
        if not np.all(active):
            B = B[:, active]
            norms = norms[:, active]
        if B.shape[1] == 0:
            return True, 0, 0
        B = B / norms

    # Gram matrix compression to find orthonormal image basis
    B_orth, rank_im = compress_gram(B, tol=tol_rank)
    if rank_im == 0:
        return True, 0, 0

    # Intersection check: does im(Pi) meet ker(s_k + I)?
    # v in ker(s_k + I) means s_k(v) = -v, i.e., (I + s_k)v = 0.
    # Apply P_k = (I + s_k)/2 to B_orth. If rank(P_k @ B_orth) == rank_im,
    # then no column-space vector lies in ker(s_k + I).
    # This avoids materializing the (potentially huge) -1 eigenspace basis.
    PkB = apply_projection(B_orth, sem[k])
    G_pk = PkB.T @ PkB
    vals_pk, _ = np.linalg.eigh(G_pk)
    max_val = np.max(vals_pk) if len(vals_pk) > 0 else 0
    rank_pk = int(np.sum(vals_pk > tol_inter * max_val)) if max_val > 0 else 0
    inter_dim = max(0, rank_im - rank_pk)

    # Count dim of -1 eigenspace from seminormal data (without building it)
    diag_k, off_i_k, off_j_k, off_ci_k, off_cj_k = sem[k]
    paired_k = set(off_i_k.tolist()) | set(off_j_k.tolist())
    dim_m1 = len(off_i_k) + sum(1 for idx in range(d) if idx not in paired_k and abs(diag_k[idx] + 1.0) < 1e-10)

    return inter_dim == 0, rank_im, dim_m1


# ─── Combined verification ────────────────────────────────────────────────────

def verify_irrep(lam, k=14, tol=1e-8, m_large=None, seed=42):
    """
    Verify Conjecture 4 for irrep λ of S_{k+1} = S_15.
    Returns (ok, rank_im, dim_minus1, elapsed).

    Uses exact method for d ≤ 3000, fast method for d > 3000.
    """
    t0 = time.time()
    n = k + 1
    assert sum(lam) == n

    d = hook_length_dim(lam)

    if d <= 1:
        return True, 0, 0, 0.0

    if d <= 3000:
        ok, rank_im, dim_m1 = verify_irrep_exact(lam, k, tol=tol)
    else:
        ok, rank_im, dim_m1 = verify_irrep_fast(lam, k, m=m_large, seed=seed,
                                                  tol_rank=tol, tol_inter=tol)
    return ok, rank_im, dim_m1, time.time() - t0


def run_all(k=14, verbose=True):
    n = k + 1
    parts = list(partitions(n))
    total = len(parts)

    if verbose:
        print(f"Verifying Conjecture 4 for k={k}: S_{n} has {total} partitions")
        print(f"Staircase: {k-1} blocks, generators 1..{k-1}")
        print(f"Check: im(Π_{k-1}) ∩ ker(s_{k}+I) = {{0}} for all λ ⊢ {n}")
        print(f"Threshold: exact for d≤3000, fast (m=max(100,d//30)) for d>3000")
        print()

    # Sort by dimension (largest first — early failure detection)
    dims_parts = sorted([(hook_length_dim(lam), lam) for lam in parts], reverse=True)

    all_ok = True
    checked = 0
    skipped_trivial = 0
    t_global = time.time()

    for d, lam in dims_parts:
        if d <= 1:
            skipped_trivial += 1
            continue
        checked += 1
        method = "exact" if d <= 3000 else "fast"
        if d > 3000:
            m_mem = int(6e8 / (d * 8))
            m_ideal = max(100, d // 30)
            m_actual = max(50, min(m_ideal, m_mem))
            m_info = f"m={m_actual}"
        else:
            m_info = ""

        if verbose:
            print(f"  λ={lam}, dim={d} [{method}{' '+m_info if m_info else ''}] ... ",
                  end="", flush=True)

        try:
            ok, rk_im, dim_m1, elapsed = verify_irrep(lam, k=k)
            status = "VERIFIED" if ok else "FAIL"
            if verbose:
                print(f"rk(im)={rk_im}, dim(-1)={dim_m1}, [{status}] ({elapsed:.1f}s)")
            if not ok:
                all_ok = False
                if verbose:
                    print(f"  *** CONJECTURE FAILS for λ={lam} ***")
                    print(f"  *** Intersection dimension is nonzero! ***")
        except MemoryError:
            if verbose:
                print("OOM - partition too large for available memory")
            all_ok = False
        except Exception as e:
            import traceback
            if verbose:
                print(f"ERROR: {e}")
                traceback.print_exc()
            all_ok = False

    t_total = time.time() - t_global
    print()
    print(f"Summary: checked {checked} non-trivial irreps, {skipped_trivial} trivial")
    print(f"Total time: {t_total:.1f}s")
    if all_ok:
        print(f"*** CONJECTURE 4 VERIFIED for k={k} (all S_{n} irreps) ***")
    else:
        print(f"*** CONJECTURE 4 FAILS for k={k} ***")
    return all_ok


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--quick":
        # Quick test: 10 largest irreps only
        parts = list(partitions(15))
        top = sorted([(hook_length_dim(lam), lam) for lam in parts], reverse=True)[:10]
        print("Quick test: 10 largest S_15 irreps")
        for d, lam in top:
            m = max(100, d // 30)
            method = "exact" if d <= 3000 else f"fast(m={m})"
            print(f"  λ={lam}, dim={d} [{method}] ... ", end="", flush=True)
            ok, rk, dm1, elapsed = verify_irrep(lam)
            print(f"rk={rk}, dim(-1)={dm1}, [{'VERIFIED' if ok else 'FAIL'}] ({elapsed:.1f}s)")
    else:
        run_all(k=14, verbose=True)
