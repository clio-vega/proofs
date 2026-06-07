"""
TASK 1b -- the CORRECT single-prime local witness for (diamond): no root mod p.

WHY this replaces the Newton-polygon hunt (task1_newton.py):
  Every Q~_b is MONIC (leading coeff 1) and primitive, so for EVERY prime p the
  p-adic Newton polygon touches height 0 at the top vertex (deg, 0) and at some
  interior unit coefficient -> it always has a slope-0 (integer) edge.  A slope-0
  edge means "roots that are p-adic units", which can perfectly well be rational.
  Hence the Newton polygon is structurally INCAPABLE of certifying "no rational
  root" for a monic polynomial.  (Confirmed: 0/110 (b,p) pairs gave a witness.)

  The right tool: since Q~_b is MONIC with integer coefficients,
        rational root  <=>  integer root  (rational root theorem, denom | LC = 1),
        and an integer root reduces to a root in F_p for every prime p.
  Therefore:
        if  Q~_b  has NO root in F_p  for SOME prime p   ==>   Q~_b has no rational root.
  This is a genuine single-prime local obstruction, and the right object to make
  congruence-uniform across the family.

This script:
  A. DIRECT complete proof of (diamond) for an extended range b=5..40 via the
     rational root theorem (test integer divisors of the constant term).  Extends
     the b<=24 record.
  B. For each b, the set of primes p<=PMAX for which Q~_b has NO root mod p
     (a local certificate).  Smallest such prime per b.
  C. Hunt for a congruence-uniform witness prime: a prime p that has no root mod p
     for ALL b in a residue class mod 4 / 8 / 12.  This is the datum (diamond) wants.
"""
import sympy as sp
from sympy import symbols, Poly, factorint, divisors, primerange, GF
import csv, os
from collections import defaultdict
from task1_newton import factor_structure   # reuse the verified Q~_b reconstruction

m = symbols('m', real=True)
RESULTS = "/home/clio/projects/code/results"

def rational_roots_monic(intpoly):
    """All rational (=integer, since monic) roots of a monic integer Poly."""
    coeffs = intpoly.all_coeffs()
    assert coeffs[0] == 1 or coeffs[0] == -1, ("not monic", coeffs[0])
    ct = int(coeffs[-1])
    if ct == 0:
        # 0 is a root; factor it out for the rest (won't happen for our Q~)
        cand = [0]
    else:
        ds = divisors(abs(ct))
        cand = sorted(set([d for d in ds] + [-d for d in ds]))
    roots = [r for r in cand if intpoly.eval(r) == 0]
    return roots

def has_root_mod_p(intpoly, p):
    """True if Q~_b (mod p) has a root in F_p (brute force over residues)."""
    coeffs = [int(c) % p for c in intpoly.all_coeffs()]
    deg = len(coeffs) - 1
    for x in range(p):
        val = 0
        for c in coeffs:
            val = (val * x + c) % p
        if val == 0:
            return True
    return False

def main():
    BMAX = 40
    PMAX = 60
    print("=" * 78)
    print(f"PART A -- DIRECT proof of (diamond): rational roots of Q~_b, b=5..{BMAX}")
    print("         (monic => rational root = integer divisor of constant term)")
    print("=" * 78)
    Qtildes = {}
    any_root = False
    for b in range(5, BMAX + 1):
        info = factor_structure(b)
        Qt = info['Qtilde']
        Qtildes[b] = Qt
        rr = rational_roots_monic(Qt)
        flag = "  <-- HAS RATIONAL ROOT" if rr else ""
        if rr:
            any_root = True
        # only print a compact line; full detail in CSV
        if b <= 16 or rr:
            print(f"  b={b:2d}  deg Q~={Qt.degree():2d}  rational roots: {rr}{flag}")
    print(f"\n  VERDICT PART A: (diamond) [Q~_b has no rational root] holds for ALL "
          f"b=5..{BMAX}: {not any_root}")
    print(f"  (extends the b<=24 record to b<={BMAX})")

    print("\n" + "=" * 78)
    print(f"PART B -- single-prime local witnesses: primes p<= {PMAX} with NO root mod p")
    print("=" * 78)
    no_root_primes = {}      # b -> sorted list of witness primes
    for b in range(5, BMAX + 1):
        Qt = Qtildes[b]
        wit = [p for p in primerange(2, PMAX + 1) if not has_root_mod_p(Qt, p)]
        no_root_primes[b] = wit
        if b <= 16:
            sm = wit[0] if wit else None
            print(f"  b={b:2d} (b%4={b%4})  smallest no-root prime: {sm}   "
                  f"all <= {PMAX}: {wit[:12]}")

    print("\n" + "=" * 78)
    print("PART C -- congruence-uniform witness hunt")
    print("=" * 78)
    # which b's does each prime witness?
    prime_hits = defaultdict(set)
    for b in range(5, BMAX + 1):
        for p in no_root_primes[b]:
            prime_hits[p].add(b)
    allb = set(range(5, BMAX + 1))
    uniform = [p for p in prime_hits if prime_hits[p] == allb]
    print(f"  prime with NO root mod p for ALL b=5..{BMAX}: {uniform if uniform else 'NONE'}")
    for mod in (4, 8, 12):
        print(f"\n  -- by residue b mod {mod} --")
        classes = defaultdict(set)
        for b in range(5, BMAX + 1):
            classes[b % mod].add(b)
        for res in sorted(classes):
            bs = classes[res]
            good = sorted(p for p in prime_hits if bs <= prime_hits[p])
            print(f"     b%{mod}={res}: class-uniform no-root primes = {good[:10] if good else 'NONE'}")

    # CSV
    csv_path = os.path.join(RESULTS, "Qb-modp-witnesses.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["b", "b_mod4", "deg_Qtilde", "rational_roots",
                    "smallest_no_root_prime", "no_root_primes_le_%d" % PMAX])
        for b in range(5, BMAX + 1):
            Qt = Qtildes[b]
            rr = rational_roots_monic(Qt)
            wit = no_root_primes[b]
            w.writerow([b, b % 4, Qt.degree(), rr,
                        wit[0] if wit else "", ";".join(map(str, wit))])
    print(f"\n  wrote {csv_path}")

    # cross-check: every "no root mod p" prime must indeed give no integer root
    print("\n  [cross-check] for each b, confirm a no-root-mod-p prime is consistent")
    ok = True
    for b in range(5, BMAX + 1):
        rr = rational_roots_monic(Qtildes[b])
        for p in no_root_primes[b]:
            for r in rr:
                if Qtildes[b].eval(r) % p == 0:   # would contradict no-root-mod-p
                    ok = False
    print(f"  consistency (no integer root reduces into a no-root-mod-p prime): {ok}")
    return Qtildes, no_root_primes, prime_hits

if __name__ == "__main__":
    main()
