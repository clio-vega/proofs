"""
TASK 1c -- finer-congruence witness hunt for the HARD classes b = 2,3 mod 4.

task1b established:
  * (diamond) holds for all b=5..40 (no rational root);
  * p=2 is a class-uniform "no root mod p" witness for b = 0,1 mod 4 (the 2-adic
    classes already PROVED in memory);
  * for b = 2,3 mod 4, NO single prime <=60 is uniform over the whole mod-4 class,
    BUT b = 10 mod 12 admitted p=47.

That hint says: split the hard classes by a FINER modulus M and look for a
congruence-class-dependent witness prime p = p(b mod M).  This script:
  - extends to b=5..BMAX, primes p<=PMAX;
  - for moduli M in {4,8,12,16,24}, reports for EACH residue class (with >=MINMEMB
    members) the set of primes that have no root mod p for every b in that class;
  - flags whether every hard class (b=2,3 mod4) is covered by SOME M.
"""
import sympy as sp
from sympy import primerange
from collections import defaultdict
import csv, os
from task1_newton import factor_structure
from task1b_modp import has_root_mod_p, rational_roots_monic

RESULTS = "/home/clio/projects/code/results"
BMAX = 64
PMAX = 120
MINMEMB = 4     # require a residue class to have at least this many members to trust it

def main():
    print(f"Building Q~_b for b=5..{BMAX} ...", flush=True)
    Qt = {}
    for b in range(5, BMAX + 1):
        Qt[b] = factor_structure(b)['Qtilde']
        if b % 8 == 0:
            print(f"   ... b={b}", flush=True)

    # (diamond) for b<=40 already proved rigorously in task1b (rational-root theorem);
    # here we focus purely on the congruence-uniform mod-p witness hunt, which does
    # NOT require factoring the (huge) constant terms.

    # no-root-mod-p witness sets
    print(f"computing no-root-mod-p witnesses (p<={PMAX}) ...", flush=True)
    no_root = {}
    for b in range(5, BMAX + 1):
        no_root[b] = set(p for p in primerange(2, PMAX + 1) if not has_root_mod_p(Qt[b], p))

    rows = []
    for mod in (4, 8, 12, 16, 24):
        print("\n" + "=" * 78)
        print(f"  modulus M = {mod}")
        print("=" * 78)
        classes = defaultdict(list)
        for b in range(5, BMAX + 1):
            classes[b % mod].append(b)
        for res in sorted(classes):
            bs = classes[res]
            if len(bs) < MINMEMB:
                continue
            common = set.intersection(*[no_root[b] for b in bs])
            smallest = min(common) if common else None
            hard = (res % 4) in (2, 3)
            tag = "HARD" if hard else "easy"
            rows.append(dict(mod=mod, res=res, n_members=len(bs),
                             members=bs, witness_primes=sorted(common),
                             smallest=smallest, hard=hard))
            print(f"   b%{mod}={res:2d} [{tag}] ({len(bs)} members {bs[:8]}{'...' if len(bs)>8 else ''}): "
                  f"uniform witness primes = {sorted(common)[:8] if common else 'NONE'}")

    # ---- the key question: is every hard class covered at SOME modulus? ----
    print("\n" + "=" * 78)
    print("  COVERAGE of the hard classes b = 2,3 mod 4")
    print("=" * 78)
    # For each fine modulus, list hard sub-classes and whether covered
    for mod in (4, 8, 12, 16, 24):
        covered, uncovered = [], []
        classes = defaultdict(list)
        for b in range(5, BMAX + 1):
            classes[b % mod].append(b)
        for res in sorted(classes):
            if (res % 4) not in (2, 3):
                continue
            bs = classes[res]
            if len(bs) < MINMEMB:
                continue
            common = set.intersection(*[no_root[b] for b in bs])
            (covered if common else uncovered).append((res, sorted(common)[:3]))
        print(f"  M={mod:2d}: covered hard residues {[(r,p) for r,p in covered]}")
        if uncovered:
            print(f"         UNCOVERED hard residues {[r for r,_ in uncovered]}")

    # CSV
    csv_path = os.path.join(RESULTS, "Qb-congruence-witnesses.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["mod", "res", "n_members", "members",
                                           "witness_primes", "smallest", "hard"])
        w.writeheader()
        for r in rows:
            rr = dict(r)
            rr["members"] = ";".join(map(str, r["members"]))
            rr["witness_primes"] = ";".join(map(str, r["witness_primes"]))
            w.writerow(rr)
    print(f"\n  wrote {csv_path}")

if __name__ == "__main__":
    main()
