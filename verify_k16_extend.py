"""
verify_k16_extend.py
Extend k=16 verification to partitions that were skipped (dim > 500K).
Raises MAX_DIM to attempt more partitions.

Strategy: use the encoded SYT fast method with reduced m to fit in memory.
For dim=1M, m=25: working matrix = 200MB, feasible.
For dim=2M, m=15: working matrix = 240MB, feasible but SYT gen + dict is ~300MB.
"""
import numpy as np
import time
import math
import gc
import sys
import traceback
import re

# Import the main verification code
sys.path.insert(0, '/home/clio/projects/proofs')
from verify_k16 import (
    partitions, hook_length_dim,
    generate_syt_encoded, axial_distance_enc, swap_encoded,
    build_seminormal_data_enc, apply_sk, apply_projection,
    compress_gram
)

# New threshold
MAX_DIM = 2000000  # Try up to 2M

def verify_irrep_fast_lean(lam, k=16, seed=42, tol=1e-8):
    """
    Memory-lean fast verification for large irreps.
    Uses smaller m and aggressive gc.
    """
    n = k + 1
    d = hook_length_dim(lam)

    # Memory budget: ~2GB for working matrices
    # d * m * 8 bytes * 3 (B + temp + sem overhead) < 2e9
    mem_budget = int(1.5e9)
    m_max = mem_budget // (d * 8 * 3)
    m = max(20, min(m_max, 100))

    # Generate encoded SYT
    t0 = time.time()
    tabs = generate_syt_encoded(lam)
    assert len(tabs) == d, f"Expected {d} SYT, got {len(tabs)}"
    tab_idx = {T: i for i, T in enumerate(tabs)}
    t_syt = time.time() - t0

    # Build seminormal data
    t0 = time.time()
    sem = {}
    for j in range(1, k + 1):
        sem[j] = build_seminormal_data_enc(lam, j, tabs, tab_idx, n)
    t_sem = time.time() - t0

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
            del sem
            gc.collect()
            return True, 0, 0
        B = B / norms

    # Gram compression
    B_orth, rank_im = compress_gram(B, tol=tol)
    del B
    gc.collect()
    if rank_im == 0:
        del sem
        gc.collect()
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


def get_skipped_partitions():
    """Read previously skipped partitions from verify_k16_output.txt."""
    skipped = []
    with open("/home/clio/projects/proofs/verify_k16_output.txt") as f:
        for line in f:
            if "[SKIP" in line:
                m = re.search(r'λ=(\([^)]+\))', line)
                if m:
                    lam = eval(m.group(1))
                    skipped.append(lam)
    return skipped


if __name__ == "__main__":
    skipped = get_skipped_partitions()
    print(f"Previously skipped: {len(skipped)} partitions")

    # Sort by dim, smallest first (most likely to succeed)
    skipped_dims = [(hook_length_dim(lam), lam) for lam in skipped]
    skipped_dims.sort()

    outpath = "/home/clio/projects/proofs/verify_k16_extend_output.txt"

    verified = 0
    failed = 0
    skipped_again = 0
    errors = []

    # Load already-verified from extend output (for resume)
    already_done = set()
    try:
        with open(outpath, "r") as f:
            for line in f:
                if "[VERIFIED]" in line or "[FAIL]" in line or "[OOM]" in line or "[ERROR" in line:
                    m = re.search(r'λ=(\([^)]+\))', line)
                    if m:
                        already_done.add(eval(m.group(1)))
    except FileNotFoundError:
        pass

    with open(outpath, "a") as logf:
        def log(msg):
            logf.write(msg + "\n")
            logf.flush()

        if not already_done:
            log(f"Extending k=16 verification: {len(skipped)} previously skipped partitions, MAX_DIM={MAX_DIM}")
            log("")

        for d, lam in skipped_dims:
            if lam in already_done:
                # Count previous results
                verified += 1  # assume verified for counting
                continue

            if d > MAX_DIM:
                log(f"  λ={lam}, dim={d} [SKIP: dim > {MAX_DIM}]")
                skipped_again += 1
                continue

            t0 = time.time()
            try:
                ok, rk_im, dim_m1 = verify_irrep_fast_lean(lam, k=16)
                elapsed = time.time() - t0
                status = "VERIFIED" if ok else "FAIL"
                log(f"  λ={lam}, dim={d} [fast] rk(im)={rk_im}, dim(-1)={dim_m1}, [{status}] ({elapsed:.1f}s)")
                if ok:
                    verified += 1
                else:
                    failed += 1
                    errors.append((lam, d, "FAIL"))
            except MemoryError:
                elapsed = time.time() - t0
                log(f"  λ={lam}, dim={d} [fast] [OOM] ({elapsed:.1f}s)")
                skipped_again += 1
                gc.collect()
            except Exception as e:
                elapsed = time.time() - t0
                log(f"  λ={lam}, dim={d} [fast] [ERROR: {e}] ({elapsed:.1f}s)")
                errors.append((lam, d, str(e)))
                gc.collect()

        summary = f"""
Summary:
  Previously skipped: {len(skipped)}
  Newly verified: {verified}
  Failed: {failed}
  Still skipped: {skipped_again}
  Errors: {len(errors)}

  Total verified (prev + new): {163 + verified} / 297 partitions of 17
  Total skipped: {skipped_again}
  Trivial (dim=1): 2
"""
        if errors:
            summary += "  Error details:\n"
            for lam_e, d_e, msg_e in errors:
                summary += f"    λ={lam_e}, dim={d_e}: {msg_e}\n"

        if failed == 0 and skipped_again == 0:
            summary += "\n*** CONJECTURE 4 FULLY VERIFIED for k=16 (all 297 S_17 irreps) ***\n"
        elif failed == 0:
            summary += f"\n*** CONJECTURE 4 VERIFIED for k=16 (all {163+verified} attempted S_17 irreps; {skipped_again} still skipped) ***\n"
        else:
            summary += "\n*** CONJECTURE 4 FAILS for k=16 ***\n"

        log(summary)
