"""JOB 3 / JOB 4 finite check — no integer root of Q_b in [b, floor(0.33 b^2)].

Q_b is monic (verified), so every rational root is an integer dividing Q_b(0).
We don't factor the (huge) constant term; instead we modular-sieve the window.

For each of several primes p, reduce Q_b mod p and find all residues r with
Q_b(r) == 0 (mod p).  A genuine integer root m in the window must satisfy
m mod p in that set for EVERY p, so we intersect candidate residue classes.
Survivors (almost always none) are then checked exactly with Q_b.eval(m).

This is the finite check the no-integer-root proof leans on.  Frontier: b<=120.
"""
import sys, time, csv
import sympy as sp
from sympy import nextprime
from _qbcore import load_cache, m

PRIMES = [1000003, 1000033, 1000037, 1000039]  # a few primes ~1e6

def zero_residues(coeffs_int, p):
    """Residues r in [0,p) with poly(r)==0 mod p, via Horner over all r.
    coeffs_int: high-to-low integer coeff list."""
    cs = [c % p for c in coeffs_int]
    zeros = set()
    for r in range(p):
        v = 0
        for c in cs:
            v = (v * r + c) % p
        if v == 0:
            zeros.add(r)
    return zeros

def main():
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    cache = load_cache()
    rows = []
    print(f"{'b':>4} {'deg':>4} {'window_hi':>10} {'exact_survivors':>16} {'time':>7}")
    for b in range(2, bmax+1):
        if b % 4 not in (2, 3):
            continue
        if b not in cache:
            continue
        t0 = time.time()
        Q = cache[b]
        deg = Q.degree()
        coeffs = [int(c) for c in Q.all_coeffs()]
        lo, hi = b, (33*b*b)//100
        # Build per-prime zero residue sets, but evaluating over full [0,p) is too
        # slow for p~1e6.  Instead: window is small (<= ~4750), just sieve the window
        # directly mod a few primes -- a candidate must vanish mod ALL primes.
        candidates = None
        for p in PRIMES:
            cs = [c % p for c in coeffs]
            good = set()
            for x in range(lo, hi+1):
                xr = x % p
                v = 0
                for c in cs:
                    v = (v * xr + c) % p
                if v == 0:
                    good.add(x)
            candidates = good if candidates is None else (candidates & good)
            if not candidates:
                break
        # exact check of survivors
        exact = [x for x in sorted(candidates) if Q.eval(x) == 0]
        ok = "NONE" if not exact else f"ROOTS!{exact}"
        dt = time.time()-t0
        print(f"{b:>4} {deg:>4} {hi:>10} {ok:>16} {dt:>6.2f}s")
        if exact:
            print("  *** INTEGER ROOT IN WINDOW -- LAW WOULD FAIL ***", b, exact)
        rows.append({"b": b, "deg": deg, "window_lo": lo, "window_hi": hi,
                     "integer_roots_in_window": ";".join(map(str, exact)) or "NONE",
                     "time_s": round(dt, 2)})
        sys.stdout.flush()
    with open("results/job3_window.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    bad = [r for r in rows if r["integer_roots_in_window"] != "NONE"]
    print(f"\nWrote results/job3_window.csv ({len(rows)} rows). "
          f"{'ALL CLEAR' if not bad else 'VIOLATIONS: '+str(bad)}")

if __name__ == "__main__":
    main()
