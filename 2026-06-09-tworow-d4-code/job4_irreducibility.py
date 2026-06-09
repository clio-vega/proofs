"""JOB 4 — irreducibility census of Q_b over QQ, b = 2,3 mod 4, via mod-p factoring.

Factoring degree-~60 integer polys over QQ (Zassenhaus) is too slow.  Use the
standard fast certificate instead:

  (CERT-IRR)  If Q_b mod p is irreducible (single factor of full degree) for some
              prime p not dividing lc(Q_b)=1, then Q_b is irreducible over QQ.

For an S_d / large-Galois family a prime giving a d-cycle (irreducible reduction)
exists with positive density (Frobenius), so a small p usually works.  When no
single p gives full irreducibility we fall back to the DEGREE-COMPATIBILITY test:
any QQ-factor degree must be a subset-sum of the mod-p irreducible-factor degrees
for EVERY good p; intersect over many p.  If only {0, deg} survive -> irreducible.

Also: is disc(Q_b) a perfect square? (Galois in A_d iff yes.)

Writes results/job4_irred.csv incrementally.
"""
import sys, time, csv
import sympy as sp
from sympy import integer_nthroot, primerange
from _qbcore import load_cache, m

def is_square(n):
    n = int(n)
    if n < 0:
        return False
    _, exact = integer_nthroot(n, 2)
    return exact

def modp_factor_degrees(coeffs, p):
    """Irreducible factor degrees of Q_b mod p (with multiplicity)."""
    P = sp.Poly(coeffs, m, modulus=p)
    if P.degree() < len(coeffs)-1:   # leading coeff vanished mod p; skip
        return None
    fl = P.factor_list()
    degs = []
    for f, mult in fl[1]:
        degs += [f.degree()] * mult
    return sorted(degs)

def subset_sums(degs):
    sums = {0}
    for d in degs:
        sums |= {s + d for s in sums}
    return sums

def main():
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    cache = load_cache()
    primes = list(primerange(3, 3000))
    rows = []
    exceptions = []
    csvf = open("results/job4_irred.csv", "w", newline="")
    writer = None
    print(f"{'b':>4} {'deg':>4} {'irred?':>8} {'how':>22} {'disc_sq?':>9} {'t':>6}")
    for b in sorted(cache):
        if b > bmax:
            continue
        Q = cache[b]
        coeffs = [int(c) for c in Q.all_coeffs()]
        deg = Q.degree()
        t0 = time.time()
        irred = None
        how = ""
        # try single-prime full irreducibility
        compat = None
        witness_p = None
        for p in primes:
            degs = modp_factor_degrees(coeffs, p)
            if degs is None:
                continue
            if len(degs) == 1 and degs[0] == deg:
                irred = True
                how = f"irred mod {p}"
                witness_p = p
                break
            ss = subset_sums(degs)
            achievable = {s for s in ss if 0 < s < deg}
            compat = achievable if compat is None else (compat & achievable)
            if compat is not None and not compat:
                irred = True
                how = "deg-incompat (multi-p)"
                break
        if irred is None:
            # exhausted primes without certifying; record remaining compat set
            irred = False
            how = f"UNCERT compat={sorted(compat)[:8] if compat else compat}"
        disc = int(sp.discriminant(Q.as_expr(), m))
        sq = is_square(disc)
        dt = time.time()-t0
        line = (f"{b:>4} {deg:>4} {('YES' if irred else 'NO?'):>8} {how:>22} "
                f"{('yes' if sq else 'no'):>9} {dt:>5.1f}s")
        print(line, flush=True)
        if not irred:
            exceptions.append(b)
            print(f"   *** NOT CERTIFIED IRREDUCIBLE b={b} -- investigate ***")
        row = {"b": b, "deg": deg, "irreducible_certified": irred,
               "method": how, "disc_is_square": sq, "time_s": round(dt, 2)}
        rows.append(row)
        if writer is None:
            writer = csv.DictWriter(csvf, fieldnames=list(row.keys()))
            writer.writeheader()
        writer.writerow(row); csvf.flush()
    csvf.close()
    print(f"\n{len(rows)} rows -> results/job4_irred.csv")
    nirr = sum(1 for r in rows if r["irreducible_certified"])
    print(f"irreducible certified: {nirr}/{len(rows)}")
    if exceptions:
        print("UNCERTIFIED b:", exceptions)
    else:
        print("ALL Q_b CERTIFIED IRREDUCIBLE over QQ, b<=", bmax)
    nsq = [r["b"] for r in rows if r["disc_is_square"]]
    print(f"disc square at b in {nsq if nsq else '[] (none -> Galois never in A_d)'}")

if __name__ == "__main__":
    main()
