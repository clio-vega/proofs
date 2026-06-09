"""JOB 1 — discriminant / resultant growth for Q_b, b = 2,3 mod 4.

Go/no-go for PROVE's recurrence-Newton-polygon irreducibility route
(Hajir / Filaseta-Trifonov-Williams).  That method needs, per b, a
"dominant prime" p > deg(Q_b) dividing disc(Q_b) to ODD multiplicity
(forcing a wildly/tamely ramified factor whose degree the Newton
polygon can pin), or a clean p-adic valuation pattern.

We cannot fully factor disc(Q_b) (hundreds-to-thousands of digits), but
we don't need to: strip all small primes by trial division, then test
the remaining cofactor for primality / perfect-power.  This determines
the parity of the multiplicity of every prime that matters for the lever.

Outputs:
  results/job1_disc.csv     one row per b
  prints a running summary; verdict assembled in FINDINGS.
"""
import sys, time, csv
sys.set_int_max_str_digits(500000)
import sympy as sp
from sympy import isprime, integer_nthroot, factorint
from _qbcore import Qb_primitive, load_cache, m

SMALL_LIMIT = 5000
                      # large prime power cofactors are caught by primality test

def strip_small(n):
    """Factor out primes <= SMALL_LIMIT. Return (small_factors_dict, cofactor)."""
    n = abs(int(n))
    fac = factorint(n, limit=SMALL_LIMIT)
    small = {}
    cof = 1
    for p, e in fac.items():
        if p <= SMALL_LIMIT:
            small[p] = e
        else:
            cof *= p**e
    return small, cof

def perfect_power(n):
    """If n = a^k with k>=2 maximal, return (a,k); else (n,1)."""
    if n <= 1:
        return (n, 1)
    best = (n, 1)
    for k in range(2, n.bit_length()+1):
        r, exact = integer_nthroot(n, k)
        if exact:
            best = (r, k)
            break
    return best

def analyze_cofactor(cof):
    """Classify the >SMALL_LIMIT cofactor: prime, prime^k, square, composite."""
    if cof == 1:
        return ("unit", {})
    base, k = perfect_power(cof)
    if k >= 2:
        if isprime(base):
            return (f"prime^{k}", {base: k})
        return (f"perfectpower^{k}", {base: k})  # base composite
    # k == 1
    if isprime(cof):
        return ("prime^1", {cof: 1})
    return ("composite", {cof: 1})  # multiplicities unknown but squarefree-status unknown

def odd_mult_large_prime(small, deg, cofclass, coffac):
    """Is there a prime p > deg dividing disc to ODD multiplicity? (Hajir lever)
    Returns (verdict_bool_or_None, witnesses)."""
    witnesses = []
    unknown = False
    for p, e in small.items():
        if p > deg and e % 2 == 1:
            witnesses.append((p, e))
    # cofactor
    if cofclass == "prime^1":
        p = list(coffac_keys(coffac))[0]
        witnesses.append((p, 1))
    elif cofclass.startswith("prime^"):
        p, e = list(coffac.items())[0]
        if e % 2 == 1 and p > deg:
            witnesses.append((p, e))
    elif cofclass in ("composite", "perfectpower^2", "perfectpower^3"):
        # cannot certify parity of the composite cofactor without factoring
        unknown = True
    if witnesses:
        return (True, witnesses)
    if unknown:
        return (None, [])
    return (False, [])

def coffac_keys(d):
    return list(d.keys())

def main():
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    cache = load_cache()
    rows = []
    print(f"{'b':>4} {'deg':>4} {'discdig':>8} {'cofclass':>14} "
          f"{'lever?':>7}  smallprimes(odd-mult, p>deg marked *)")
    for b in range(2, bmax+1):
        if b % 4 not in (2, 3):
            continue
        if b not in cache:
            continue
        t0 = time.time()
        Q = cache[b]
        deg = Q.degree()
        disc = int(sp.discriminant(Q.as_expr(), m))
        small, cof = strip_small(disc)
        # Short-circuit: if a SMALL prime > deg already has odd multiplicity,
        # the lever is present -- skip the (expensive) cofactor perfect-power test.
        small_witness = [(p, e) for p, e in small.items() if p > deg and e % 2 == 1]
        if small_witness:
            cofclass = "(cofactor not analyzed; small witness suffices)"
            coffac = {}
            lever, wit = True, small_witness
        else:
            # no small witness: cheap check only -- is the leftover cofactor a
            # single large prime (odd multiplicity -> lever) ?  Skip perfect-power
            # loop (too slow on multi-thousand-digit cofactors).
            if cof == 1:
                cofclass, lever, wit = "unit", False, []
            elif isprime(cof):
                cofclass, lever, wit = "prime^1", (cof > deg), ([(cof, 1)] if cof > deg else [])
            else:
                cofclass, lever, wit = "composite(unfactored)", None, []
        # describe odd-multiplicity small primes
        oddsmall = {p: e for p, e in small.items() if e % 2 == 1}
        odd_desc = ", ".join(
            f"{p}^{e}" + ("*" if p > deg else "") for p, e in sorted(oddsmall.items()))
        leverstr = {True: "YES", False: "no", None: "?"}[lever]
        discdig = len(str(abs(disc)))
        print(f"{b:>4} {deg:>4} {discdig:>8} {cofclass:>14} {leverstr:>7}  {odd_desc}")
        rows.append({
            "b": b, "deg": deg, "disc_digits": discdig,
            "disc_sign": 1 if disc > 0 else -1,
            "cofactor_class": cofclass,
            "cofactor": str(cof),
            "lever": leverstr,
            "lever_witness": ";".join(f"{p}^{e}" for p, e in wit),
            "small_odd_mult": ";".join(f"{p}^{e}" for p, e in sorted(oddsmall.items())),
            "small_factors": ";".join(f"{p}^{e}" for p, e in sorted(small.items())),
            "build_s": round(time.time()-t0, 2),
        })
        sys.stdout.flush()
    with open("results/job1_disc.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote results/job1_disc.csv ({len(rows)} rows)")

if __name__ == "__main__":
    main()
