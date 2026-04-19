"""
verify_k16.py
Verify Conjecture 4 for k=16: for every irrep λ ⊢ 17,
    im(Π_{15}) ∩ ker(s_{16} + I) = {0}

S_17 has 297 partitions. Largest irreps have dim ~ 5 million.

Strategy:
 - EXACT (dim ≤ 3000): SVD-based image tracking
 - FAST (3000 < dim ≤ MAX_DIM): random projection + Gram compression
 - SKIP (dim > MAX_DIM): too large for this container's memory

For the FAST method, SYT are encoded as numpy int8 arrays (row of each entry)
to minimize Python object overhead. Each SYT of shape λ ⊢ 17 needs only 17 bytes
in this encoding, vs ~1KB+ for nested Python lists/tuples.

Seminormal form: s_j|T⟩ = (1/ρ)|T⟩ + √(1-1/ρ²)|T'⟩
where ρ = axial_distance(T,j) and T' = swap_j(T).
"""
import numpy as np
import time
import math
import gc
import sys
import traceback

# Maximum dimension we attempt.
# SYT enumeration + numpy arrays: d=500K needs ~500K*17 bytes (SYT) + d*m*8 (work matrix)
# At d=500K, m=50: SYT=8.5MB, work=200MB. Feasible.
# At d=800K, m=40: SYT=13.6MB, work=256MB. Feasible for arrays, but SYT generation time matters.
# The real bottleneck is Python's recursive SYT enumeration: ~500K/sec for n=17.
# d=1M takes ~2 sec, d=5M takes ~10 sec for enumeration alone — actually manageable.
# But we need a dict for lookups: d=1M entries * ~100 bytes overhead = 100MB.
# Let's try MAX_DIM = 500000 first.
MAX_DIM = 500000


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


# ─── Standard tableaux with compact encoding ─────────────────────────────────

def generate_syt_encoded(lam):
    """
    Generate all SYT of shape lam. Each SYT is encoded as a bytes object
    of length n, where byte[v-1] = row of entry v. This is extremely compact
    and hashable for dict lookups.

    Returns: list of bytes objects, one per SYT.
    """
    n = sum(lam)
    if n == 0:
        return [b'']

    results = []
    rows = bytearray(n)  # rows[v-1] = row of entry v
    row_lens = [0] * len(lam)
    # We also need column info for axial distance, but we can reconstruct it
    # from row_lens at placement time. Store (row, col) but encode efficiently.
    # Actually, let's store row positions as bytes AND column positions as bytes.
    cols = bytearray(n)

    def fill(v):
        if v > n:
            # Encode as (rows_bytes, cols_bytes) combined into single bytes
            results.append(bytes(rows) + bytes(cols))
            return
        for i in range(len(lam)):
            if row_lens[i] < lam[i]:
                if i == 0 or row_lens[i] < row_lens[i - 1]:
                    c = row_lens[i]
                    rows[v - 1] = i
                    cols[v - 1] = c
                    row_lens[i] += 1
                    fill(v + 1)
                    row_lens[i] -= 1

    fill(1)
    return results


def axial_distance_enc(enc, k, n):
    """Axial distance for encoded SYT. k is 1-indexed."""
    r1, c1 = enc[k - 1], enc[n + k - 1]
    r2, c2 = enc[k], enc[n + k]
    return (c2 - r2) - (c1 - r1)


def swap_encoded(enc, k, n, lam):
    """
    Swap entries k and k+1 in encoded SYT. Returns new encoding if valid, else None.
    """
    enc_new = bytearray(enc)
    # Swap row/col positions of entries k and k+1
    enc_new[k - 1], enc_new[k] = enc[k], enc[k - 1]
    enc_new[n + k - 1], enc_new[n + k] = enc[n + k], enc[n + k - 1]

    # Check validity: for each value v, row[v] and col[v] must be consistent
    # with a standard tableau. We need:
    # 1. For consecutive values in the same row: left < right
    # 2. For consecutive values in the same column: top < bottom
    # Equivalently: for each v, if v-1 is at (r',c') and v is at (r,c):
    #   - if r == r', then c > c' (same row, must be right)
    #   - if c == c', then r > r' (same column, must be below)
    # But this isn't sufficient. Need to check the grid structure.

    # Build a grid and verify
    grid = {}
    for v in range(n):
        r, c = enc_new[v], enc_new[n + v]
        grid[(r, c)] = v + 1

    for v in range(n):
        r, c = enc_new[v], enc_new[n + v]
        val = v + 1
        if c > 0 and (r, c - 1) in grid and grid[(r, c - 1)] >= val:
            return None
        if r > 0 and (r - 1, c) in grid and grid[(r - 1, c)] >= val:
            return None

    return bytes(enc_new)


# ─── Seminormal form ──────────────────────────────────────────────────────────

def build_seminormal_data_enc(lam, k, tabs, tab_idx, n):
    """Build sparse seminormal representation using encoded tableaux."""
    if k >= n:
        d = len(tabs)
        return (np.zeros(d), np.array([], dtype=np.intp),
                np.array([], dtype=np.intp), np.array([]), np.array([]))
    d = len(tabs)
    diag = np.zeros(d)
    oi, oj, oci, ocj = [], [], [], []
    seen = set()

    for idx_i in range(d):
        if idx_i in seen:
            continue
        T = tabs[idx_i]
        rho = axial_distance_enc(T, k, n)

        T_sw = swap_encoded(T, k, n, lam)
        if T_sw is not None and T_sw in tab_idx:
            j = tab_idx[T_sw]
            c = math.sqrt(max(0.0, 1.0 - 1.0 / rho**2))
            oi.append(idx_i)
            oj.append(j)
            oci.append(1.0 / rho)
            ocj.append(c)
            seen.add(idx_i)
            seen.add(j)
            diag[idx_i] = 1.0 / rho
            diag[j] = -1.0 / rho
            continue
        seen.add(idx_i)
        diag[idx_i] = 1.0 / rho

    return (diag,
            np.array(oi, dtype=np.intp),
            np.array(oj, dtype=np.intp),
            np.array(oci),
            np.array(ocj))


# Standard list-based SYT (for exact method on small dims)
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

def axial_distance_std(T, k):
    res1 = find_entry(T, k)
    res2 = find_entry(T, k + 1)
    if res1 is None or res2 is None:
        return None
    return (res2[1] - res2[0]) - (res1[1] - res1[0])

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

def build_seminormal_data_std(lam, k, tabs, tab_idx):
    n = sum(lam)
    if k >= n:
        d = len(tabs)
        return (np.zeros(d), np.array([], dtype=np.intp),
                np.array([], dtype=np.intp), np.array([]), np.array([]))
    d = len(tabs)
    diag = np.zeros(d)
    oi, oj, oci, ocj = [], [], [], []
    seen = set()

    for T in tabs:
        i = tab_idx[tab_key(T)]
        if i in seen:
            continue
        rho = axial_distance_std(T, k)
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


# ─── Apply operators ─────────────────────────────────────────────────────────

def apply_sk(B, sem):
    diag, off_i, off_j, off_ci, off_cj = sem
    sB = diag[:, np.newaxis] * B
    if len(off_i) > 0:
        sB[off_i] += off_cj[:, np.newaxis] * B[off_j]
        sB[off_j] += off_cj[:, np.newaxis] * B[off_i]
    return sB

def apply_projection(B, sem):
    return (B + apply_sk(B, sem)) * 0.5


# ─── Gram matrix compression ────────────────────────────────────────────────

def compress_gram(B, tol=1e-8):
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


# ─── Exact method (d ≤ 3000) ─────────────────────────────────────────────────

def build_plus1_eigenspace(d, sem_j):
    diag, off_i, off_j, off_ci, off_cj = sem_j
    basis = []
    paired = set(off_i.tolist()) | set(off_j.tolist())
    for idx in range(len(off_i)):
        i, j = off_i[idx], off_j[idx]
        rho = 1.0 / off_ci[idx]
        c = off_cj[idx]
        v = np.zeros(d)
        v[i] = c
        v[j] = 1.0 - 1.0 / rho
        nrm = np.linalg.norm(v)
        if nrm > 1e-14:
            basis.append(v / nrm)
        else:
            v2 = np.zeros(d); v2[i] = 1.0
            basis.append(v2)
    for idx in range(d):
        if idx not in paired and diag[idx] > 0:
            v = np.zeros(d); v[idx] = 1.0
            basis.append(v)
    if not basis:
        return np.zeros((d, 0))
    return np.column_stack(basis)


def apply_projection_svd(B, sem, tol=1e-8):
    PB = apply_projection(B, sem)
    U, s, _ = np.linalg.svd(PB, full_matrices=False)
    rank = int(np.sum(s > tol))
    return U[:, :rank]


def build_minus1_basis(d, sem_k):
    diag_k, off_i, off_j, off_ci, off_cj = sem_k
    minus1_vecs = []
    paired = set(off_i.tolist()) | set(off_j.tolist())
    for idx in range(len(off_i)):
        i, j = off_i[idx], off_j[idx]
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
            v = np.zeros(d); v[idx] = 1.0
            minus1_vecs.append(v)
    if not minus1_vecs:
        return np.zeros((d, 0))
    return np.column_stack(minus1_vecs)


def verify_irrep_exact(lam, k, tol=1e-8):
    n = k + 1
    tabs = standard_tableaux(lam)
    d = len(tabs)
    tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}

    sem = {}
    for j in range(1, k + 1):
        sem[j] = build_seminormal_data_std(lam, j, tabs, tab_idx)

    B = build_plus1_eigenspace(d, sem[1])
    if B.shape[1] == 0:
        return True, 0, 0

    for block in range(2, k):
        for j in range(block, 0, -1):
            B = apply_projection_svd(B, sem[j], tol=tol)
            if B.shape[1] == 0:
                return True, 0, 0

    rank_im = B.shape[1]

    # Count dim of -1 eigenspace from seminormal data
    diag_k = sem[k][0]
    off_i_k = sem[k][1]
    paired_k = set(off_i_k.tolist()) | set(sem[k][2].tolist())
    dim_m1 = len(off_i_k) + sum(1 for idx in range(d) if idx not in paired_k and abs(diag_k[idx] + 1.0) < 1e-10)

    if dim_m1 == 0:
        return True, rank_im, 0

    # Intersection check via P_k (more numerically stable than Gram residual method)
    # v in ker(s_k+I) iff P_k(v) = 0, so rank drop of P_k @ B = intersection dim
    PkB = apply_projection(B, sem[k])
    sv_pk = np.linalg.svd(PkB, compute_uv=False)
    rank_pk = int(np.sum(sv_pk > tol))
    inter_dim = max(0, rank_im - rank_pk)

    return inter_dim == 0, rank_im, dim_m1


# ─── Fast method (d > 3000) ──────────────────────────────────────────────────

def verify_irrep_fast(lam, k, m=None, seed=42, tol=1e-8):
    """
    Fast method using encoded SYT and random projection.
    Memory-efficient: SYT encoded as 34-byte strings (for n=17).
    """
    n = k + 1
    d = hook_length_dim(lam)

    # Determine m: balance between statistical power and memory
    # Working matrix: d * m * 8 bytes
    # Keep total numpy memory under 1 GB: d * m * 8 * 2 < 1e9 (factor 2 for temp copies)
    mem_budget = int(0.5e9)  # 500 MB for working matrix (leaving room for copies + seminormal)
    m_mem = mem_budget // (d * 8)
    if m is None:
        m_ideal = max(50, d // 60)
        m = max(30, min(m_ideal, m_mem, 200))

    print(f"(m={m}) ", end="", flush=True)

    # Generate encoded SYT
    t_syt = time.time()
    tabs = generate_syt_encoded(lam)
    assert len(tabs) == d, f"Expected {d} SYT, got {len(tabs)}"
    tab_idx = {T: i for i, T in enumerate(tabs)}
    t_syt = time.time() - t_syt
    print(f"(syt {t_syt:.0f}s) ", end="", flush=True)

    # Build seminormal data for all generators
    t_sem = time.time()
    sem = {}
    for j in range(1, k + 1):
        sem[j] = build_seminormal_data_enc(lam, j, tabs, tab_idx, n)
    t_sem = time.time() - t_sem
    print(f"(sem {t_sem:.0f}s) ", end="", flush=True)

    # Free tableaux
    del tabs, tab_idx
    gc.collect()

    # Random initial basis
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((d, m)).astype(np.float64)
    norms = np.linalg.norm(B, axis=0, keepdims=True)
    B /= norms

    # Apply full staircase
    for block in range(1, k):
        for j in range(block, 0, -1):
            B = apply_projection(B, sem[j])
        norms = np.linalg.norm(B, axis=0, keepdims=True)
        active = norms[0] > 1e-30
        if not np.all(active):
            B = B[:, active]
            norms = norms[:, active]
        if B.shape[1] == 0:
            return True, 0, 0
        B = B / norms

    # Gram compression
    B_orth, rank_im = compress_gram(B, tol=tol)
    del B
    gc.collect()
    if rank_im == 0:
        return True, 0, 0

    # Intersection check via P_k
    PkB = apply_projection(B_orth, sem[k])
    G_pk = PkB.T @ PkB
    del PkB
    vals_pk, _ = np.linalg.eigh(G_pk)
    max_val = np.max(vals_pk) if len(vals_pk) > 0 else 0
    rank_pk = int(np.sum(vals_pk > tol * max_val)) if max_val > 0 else 0
    inter_dim = max(0, rank_im - rank_pk)

    # Count -1 eigenspace dimension
    diag_k = sem[k][0]
    off_i_k = sem[k][1]
    paired_k = set(off_i_k.tolist()) | set(sem[k][2].tolist())
    dim_m1 = len(off_i_k) + sum(1 for idx in range(d) if idx not in paired_k and abs(diag_k[idx] + 1.0) < 1e-10)

    del B_orth, sem
    gc.collect()

    return inter_dim == 0, rank_im, dim_m1


# ─── Combined verification ────────────────────────────────────────────────────

def verify_irrep(lam, k=16, tol=1e-8, seed=42):
    t0 = time.time()
    n = k + 1
    assert sum(lam) == n
    d = hook_length_dim(lam)

    if d <= 1:
        return True, 0, 0, 0.0

    if d <= 3000:
        ok, rank_im, dim_m1 = verify_irrep_exact(lam, k, tol=tol)
    else:
        ok, rank_im, dim_m1 = verify_irrep_fast(lam, k, seed=seed, tol=tol)

    return ok, rank_im, dim_m1, time.time() - t0


def load_verified(outpath):
    """Load already-verified partitions from previous run output."""
    verified = set()
    try:
        with open(outpath, "r") as f:
            for line in f:
                if "[VERIFIED]" in line and "λ=" in line:
                    # Extract partition from line like "  λ=(8, 4, 1, 1, 1, 1, 1), dim=..."
                    start = line.index("λ=(") + 2
                    end = line.index(")", start) + 1
                    lam_str = line[start:end]
                    try:
                        lam = eval(lam_str)
                        verified.add(lam)
                    except:
                        pass
    except FileNotFoundError:
        pass
    return verified


def run_all(k=16, verbose=True, logfile=None, resume_from=None):
    n = k + 1
    parts = list(partitions(n))
    total = len(parts)

    # Load previously verified partitions for resume
    already_verified = set()
    if resume_from:
        already_verified = load_verified(resume_from)
        if already_verified and verbose:
            print(f"Resuming: {len(already_verified)} partitions already verified", flush=True)

    def log(msg, end="\n"):
        if verbose:
            print(msg, end=end, flush=True)
        if logfile:
            logfile.write(msg + (end if end != "\n" else "\n"))
            logfile.flush()

    log(f"Verifying Conjecture 4 for k={k}: S_{n} has {total} partitions")
    log(f"Staircase: {k-1} blocks, generators 1..{k-1}")
    log(f"Check: im(Π_{k-1}) ∩ ker(s_{k}+I) = {{0}} for all λ ⊢ {n}")
    log(f"Threshold: exact d≤3000, fast d≤{MAX_DIM}, skip d>{MAX_DIM}")
    if already_verified:
        log(f"Resuming with {len(already_verified)} already verified")
    log("")

    dims_parts = sorted([(hook_length_dim(lam), lam) for lam in parts], reverse=True)

    all_ok = True
    checked = 0
    resumed = 0
    skipped_trivial = 0
    skipped_large = []
    failed = []
    t_global = time.time()

    for d, lam in dims_parts:
        if d <= 1:
            skipped_trivial += 1
            continue

        if d > MAX_DIM:
            skipped_large.append((lam, d))
            log(f"  λ={lam}, dim={d} [SKIP: dim > {MAX_DIM}]")
            continue

        if lam in already_verified:
            resumed += 1
            log(f"  λ={lam}, dim={d} [ALREADY VERIFIED]")
            checked += 1
            continue

        checked += 1
        method = "exact" if d <= 3000 else "fast"

        msg = f"  λ={lam}, dim={d} [{method}] "
        if verbose:
            print(msg, end="", flush=True)
        if logfile:
            logfile.write(msg)
            logfile.flush()

        try:
            ok, rk_im, dim_m1, elapsed = verify_irrep(lam, k=k)
            status = "VERIFIED" if ok else "FAIL"
            result = f"rk(im)={rk_im}, dim(-1)={dim_m1}, [{status}] ({elapsed:.1f}s)"
            if verbose:
                print(result, flush=True)
            if logfile:
                logfile.write(result + "\n")
                logfile.flush()
            if not ok:
                all_ok = False
                failed.append((lam, d))
                log(f"  *** CONJECTURE FAILS for λ={lam} ***")
            gc.collect()
        except MemoryError:
            msg_e = "OOM - skipped"
            if verbose:
                print(msg_e, flush=True)
            if logfile:
                logfile.write(msg_e + "\n")
                logfile.flush()
            skipped_large.append((lam, d))
            gc.collect()
        except Exception as e:
            msg_e = f"ERROR: {e}"
            if verbose:
                print(msg_e, flush=True)
                traceback.print_exc()
            if logfile:
                logfile.write(msg_e + "\n")
                logfile.flush()
            failed.append((lam, d))
            gc.collect()

    t_total = time.time() - t_global
    log("")
    log(f"Summary: checked {checked} non-trivial irreps ({resumed} resumed), {skipped_trivial} trivial, {len(skipped_large)} skipped (too large/OOM)")
    if skipped_large:
        log(f"Skipped partitions ({len(skipped_large)}):")
        for lam_s, d_s in skipped_large:
            log(f"  λ={lam_s}, dim={d_s}")
    if failed:
        log(f"FAILED partitions ({len(failed)}):")
        for lam_f, d_f in failed:
            log(f"  λ={lam_f}, dim={d_f}")
    log(f"Total time: {t_total:.1f}s")
    if all_ok and not skipped_large:
        log(f"*** CONJECTURE 4 VERIFIED for k={k} (all {total} S_{n} irreps) ***")
    elif all_ok:
        log(f"*** CONJECTURE 4 VERIFIED for k={k} (all {checked} attempted S_{n} irreps; {len(skipped_large)} skipped) ***")
    else:
        log(f"*** CONJECTURE 4 FAILS for k={k} ***")
    return all_ok


if __name__ == "__main__":
    outpath = "/home/clio/projects/proofs/verify_k16_output.txt"
    prev_outpath = "/home/clio/projects/proofs/verify_k16_prev.txt"

    # Check for previous run to resume from
    import os
    resume = None
    if os.path.exists(outpath):
        # Rename old output, use it for resume
        import shutil
        shutil.copy2(outpath, prev_outpath)
        resume = prev_outpath

    with open(outpath, "w") as f:
        run_all(k=16, verbose=True, logfile=f, resume_from=resume)
    print(f"\nOutput saved to {outpath}")
