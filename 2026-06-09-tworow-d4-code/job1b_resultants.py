"""JOB 1 (second half) — cross-resultants Res(Q_b, Q_{b'}) for nearby b.

Question for the no-shared-root / structured-coprimality angle: do distinct
members of the family Q_b ever share an integer root, and is the family
"structured-coprime"?  Since each Q_b (b>=6) is irreducible of degree >= 3
(Job 4), an integer root of Q_b is already impossible -- but we still record
the resultant prime structure to see whether the family is congruence-coupled.

For consecutive b ≡ 2,3 (mod 4) pairs we compute Res(Q_b, Q_{b'}), strip small
primes, and report the magnitude / small-prime content.  gcd over QQ is checked
(must be 1 for distinct irreducibles).
"""
import sys, csv
sys.set_int_max_str_digits(500000)
import sympy as sp
from sympy import factorint
from _qbcore import load_cache, m

SMALL_LIMIT = 50000

def small_content(n):
    n = abs(int(n))
    if n == 0:
        return "0", {}
    fac = factorint(n, limit=SMALL_LIMIT)
    small = {p: e for p, e in fac.items() if p <= SMALL_LIMIT}
    cof = 1
    for p, e in fac.items():
        if p > SMALL_LIMIT:
            cof *= p**e
    return cof, small

def main():
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    cache = load_cache()
    bs = sorted(b for b in cache if b % 4 in (2, 3) and b <= bmax and b >= 6)
    rows = []
    print(f"{'b':>4} {'b2':>4} {'gcd=1?':>7} {'res_digits':>10} small-prime content")
    for i in range(len(bs)-1):
        b, b2 = bs[i], bs[i+1]
        Q1, Q2 = cache[b].as_expr(), cache[b2].as_expr()
        g = sp.gcd(sp.Poly(Q1, m), sp.Poly(Q2, m))
        gcd1 = (g.degree() == 0)
        res = int(sp.resultant(Q1, Q2, m))
        cof, small = small_content(res)
        smalldesc = ", ".join(f"{p}^{e}" for p, e in sorted(small.items())[:12])
        print(f"{b:>4} {b2:>4} {str(gcd1):>7} {len(str(abs(res))):>10} {smalldesc}")
        rows.append({"b": b, "b2": b2, "gcd_is_1": gcd1,
                     "res_digits": len(str(abs(res))),
                     "small_primes": ";".join(f"{p}^{e}" for p, e in sorted(small.items())),
                     "large_cofactor_digits": len(str(cof)) if cof != 1 else 0})
        sys.stdout.flush()
    with open("results/job1_resultants.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\nWrote results/job1_resultants.csv ({len(rows)} rows). "
          f"gcd=1 for all pairs: {all(r['gcd_is_1'] for r in rows)}")

if __name__ == "__main__":
    main()
