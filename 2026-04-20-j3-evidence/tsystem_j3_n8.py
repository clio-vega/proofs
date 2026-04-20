"""Verify the j=3 T-system / Hirota identity at n=8:
    d_7(8)(q) * d_5(8)(q)  ?=  q^|D_8| * d_6(8)(q)^2     in Z[q]
where d_k(n) = det M_R^(k)[D_n, D_n], |D_8| = 8!/2^4 = 2520.

Strategy:
  - At each integer q_0, build the integer matrices M^(k)[D_8, D_8] for k=5,6,7
    using the fast Hecke right-action (adapted from
    /home/clio/projects/scratch/multi_point_verify.py).
  - For each large prime p in a list, reduce the integer matrix mod p and use
    python-flint nmod_mat.det to compute det mod p in <10s for 2520x2520.
  - Test the identity d_7 * d_5 == q^|D_8| * d_6^2 mod each prime.

Multi-prime mod-p evidence at one q_0 effectively certifies the integer identity
d_7(q_0) * d_5(q_0) = q_0^|D_8| * d_6(q_0)^2 (a non-trivial integer would be
detected by the mod-p check with overwhelming probability), and multi-q
evidence then certifies the polynomial identity in Z[q] up to a Schwartz-Zippel
bound.

Usage:
  python3 tsystem_j3_n8.py --n 8 --qs 2,3,5 --ks 5,6,7
  python3 tsystem_j3_n8.py --n 6 --qs 2,3,5 --ks 3,4,5    # j=2 sanity check
  python3 tsystem_j3_n8.py --n 7 --qs 2,3,5 --ks 3,4,5    # j=2 cross-check at n=7
"""

import sys
import time
import argparse
from collections import defaultdict
from itertools import permutations

try:
    import flint
    HAVE_FLINT = True
except ImportError:
    HAVE_FLINT = False


# --------------------------------------------------------------------------
# Combinatorics
# --------------------------------------------------------------------------

def descent_paired(n):
    D = []
    fh = n // 2
    for w in permutations(range(1, n+1)):
        if all(w[2*i] > w[2*i+1] for i in range(fh)):
            D.append(w)
    return D


def staircase_word_partial(k):
    """Word of i's so that prod (1+T_{w_t}) = B_1 B_2 ... B_k where
    B_j = (1+T_j)(1+T_{j-1})...(1+T_1)."""
    word = []
    for block in range(1, k+1):
        for i in range(block, 0, -1):
            word.append(i)
    return word


# --------------------------------------------------------------------------
# Hecke right-action at integer q (no mod reduction during build).
# --------------------------------------------------------------------------

def apply_one_plus_Ti_right_int(coeffs, i, q_val):
    """coeffs : dict perm -> int; returns coeffs * (1 + T_i)."""
    qm1 = q_val - 1
    via_Ti = defaultdict(int)
    for w, c in coeffs.items():
        if c == 0:
            continue
        wsi = list(w)
        wsi[i-1], wsi[i] = wsi[i], wsi[i-1]
        wsi = tuple(wsi)
        if w[i-1] < w[i]:
            via_Ti[wsi] += c
        else:
            via_Ti[wsi] += q_val * c
            via_Ti[w] += qm1 * c
    result = defaultdict(int)
    for w, c in coeffs.items():
        result[w] += c
    for w, c in via_Ti.items():
        result[w] += c
    return {w: c for w, c in result.items() if c != 0}


def compute_column_int(v, k, q_val):
    coeffs = {v: 1}
    for i in staircase_word_partial(k):
        coeffs = apply_one_plus_Ti_right_int(coeffs, i, q_val)
    return coeffs


def build_M_int_rows(n, k, q_val, D):
    """Return list of m rows, each a list of m ints. M[i][j] = coeff of T_{D[i]}
    in T_{D[j]} * Π^(k)_q evaluated at q=q_val."""
    m = len(D)
    D_to_idx = {w: i for i, w in enumerate(D)}
    M = [[0] * m for _ in range(m)]
    for j, v in enumerate(D):
        col = compute_column_int(v, k, q_val)
        for u, c in col.items():
            i = D_to_idx.get(u)
            if i is not None:
                M[i][j] = c
    return M


# --------------------------------------------------------------------------
# Determinants (using flint)
# --------------------------------------------------------------------------

def det_mod_via_flint(M_rows, p):
    Mp = flint.nmod_mat(M_rows, p)
    return int(Mp.det())


# --------------------------------------------------------------------------
# Main verification routine
# --------------------------------------------------------------------------

def verify_q_value(n, q_val, D, ks, primes, log_fn):
    """Build M^(k) over Z for each k in ks, then for each prime test
    d_{ks[2]} * d_{ks[0]} ?= q^|D_n| * d_{ks[1]}^2 mod p.
    Returns: per-prime list of bool, plus the per-k integer matrices' max bit-length.
    """
    exponent = len(D)
    matrices = {}
    build_times = {}
    maxbits = {}

    for k in ks:
        log_fn(f"  k={k}: building integer matrix at q={q_val}...")
        t0 = time.time()
        M = build_M_int_rows(n, k, q_val, D)
        t1 = time.time()
        build_times[k] = t1 - t0
        mb = max((abs(x).bit_length() for r in M for x in r), default=0)
        maxbits[k] = mb
        log_fn(f"    build {build_times[k]:6.1f}s; max entry bits = {mb}")
        matrices[k] = M

    # Now test the identity at each prime.
    results = []
    dets_per_prime = {}
    for p in primes:
        log_fn(f"  prime p (~{p.bit_length()} bits): {p}")
        d_mod = {}
        for k in ks:
            t0 = time.time()
            d = det_mod_via_flint(matrices[k], p)
            t1 = time.time()
            d_mod[k] = d
            log_fn(f"    k={k}: det = {d}  ({t1-t0:.1f}s)")
        # Test identity
        k_lo, k_mid, k_hi = ks[0], ks[1], ks[2]
        lhs = (d_mod[k_hi] * d_mod[k_lo]) % p
        rhs = (pow(q_val, exponent, p) * d_mod[k_mid] * d_mod[k_mid]) % p
        ok = (lhs == rhs)
        log_fn(f"    identity d_{k_hi}*d_{k_lo} == q^{exponent}*d_{k_mid}^2 mod p? "
               f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            log_fn(f"      lhs = {lhs}")
            log_fn(f"      rhs = {rhs}")
            log_fn(f"      diff = {(lhs - rhs) % p}")
        results.append((p, ok, lhs, rhs, dict(d_mod)))
        dets_per_prime[p] = dict(d_mod)
    return results, build_times, maxbits, dets_per_prime


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--qs", type=str, default="2,3,5",
                        help="Comma-separated integer q values to test.")
    parser.add_argument("--ks", type=str, default="5,6,7",
                        help="Comma-separated k values for which to compute "
                             "det M^(k); MUST be exactly 3 (k_lo, k_mid, k_hi).")
    parser.add_argument("--primes", type=str,
                        default="9223372036854775783,9223372036854775643,9223372036854775549",
                        help="Comma-separated primes (~63 bits each).")
    parser.add_argument("--time-budget-seconds", type=float, default=1800.0,
                        help="Per-q time budget (graceful abort).")
    parser.add_argument("--log", type=str, default=None,
                        help="Optional log file to also write to.")
    args = parser.parse_args()

    if not HAVE_FLINT:
        print("ERROR: python-flint not available; pip install --break-system-packages python-flint",
              flush=True)
        sys.exit(1)

    n = args.n
    q_vals = [int(x) for x in args.qs.split(",")]
    ks = tuple(int(x) for x in args.ks.split(","))
    primes = [int(x) for x in args.primes.split(",")]

    if len(ks) != 3:
        print("ERROR: need exactly 3 ks (k_lo, k_mid, k_hi).", flush=True)
        sys.exit(1)

    log_fp = open(args.log, 'w') if args.log else None

    def log(msg):
        print(msg, flush=True)
        if log_fp:
            log_fp.write(msg + "\n")
            log_fp.flush()

    log(f"# T-system verification, n={n}, ks={ks}")
    log(f"# Primes (count={len(primes)}):")
    for p in primes:
        log(f"#   {p}  (~{p.bit_length()} bits)")
    log(f"# q values: {q_vals}")

    t_total = time.time()
    log(f"\n## Build descent-paired set D_{n} ...")
    t0 = time.time()
    D = descent_paired(n)
    log(f"|D_{n}| = {len(D)}  (built in {time.time()-t0:.2f}s)")
    log(f"# Conjecture: d_{ks[2]} * d_{ks[0]} = q^{len(D)} * d_{ks[1]}^2  in Z[q]")

    overall = {}
    for q_val in q_vals:
        log(f"\n========================================")
        log(f"q = {q_val}")
        log(f"========================================")
        t_q0 = time.time()
        try:
            results, bt, mb, dpp = verify_q_value(n, q_val, D, ks, primes, log)
        except Exception as e:
            import traceback
            log(f"  ERROR at q={q_val}: {e}")
            log(traceback.format_exc())
            overall[q_val] = "error"
            continue
        elapsed = time.time() - t_q0
        all_pass = all(r[1] for r in results)
        log(f"  Total time at q={q_val}: {elapsed:.1f}s")
        log(f"  Per-prime results: {[r[1] for r in results]}")
        log(f"  Overall at q={q_val}: {'PASS' if all_pass else 'FAIL'}")
        overall[q_val] = "pass" if all_pass else "fail"

        if elapsed > args.time_budget_seconds:
            log(f"\n# Time budget {args.time_budget_seconds}s exceeded at q={q_val}; aborting further qs.")
            break

    log(f"\n# ============================== SUMMARY ==============================")
    log(f"# n={n}, |D_n|={len(D)}")
    log(f"# Conjecture: d_{ks[2]} * d_{ks[0]} = q^{len(D)} * d_{ks[1]}^2")
    for q_val, status in overall.items():
        log(f"#   q={q_val}: {status.upper()}")
    log(f"# Total elapsed: {time.time()-t_total:.1f}s")

    if log_fp:
        log_fp.close()


if __name__ == "__main__":
    main()
