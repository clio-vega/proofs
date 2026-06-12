"""
Job A.1 -- full table of (lambda, {val(j)}, J*, |J*|, V) for all lambda |- 2m, m <= MMAX.

Checks:
  * |J*| even whenever |J*| >= 2  (extend the m<=7 empirical wall)
  * record any |J*| = 8 cases; check whether the J* offsets stay a 2-adic box
  * **GUARD**: if a single ODD |J*| >= 3 appears, FLAG IT LOUDLY -- it breaks the
    leading-pi dichotomy and is the headline.
  * record parity of J* members (all same parity expected -> tells us a(s) vs b(s)).
"""
import sys, json
from collections import Counter
from job_jstar_engine import analyze, partitions, v2

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16

def offsets_are_2adic_box(Jstar):
    """J* = j0 + S, S = subset-sums of distinct powers of 2 (affine 2-adic box).  Since the
    generators are distinct powers of 2, S is an XOR-closed subspace; detect via XOR closure and
    extract a GF(2) basis (the generators).  Returns (is_box, generators)."""
    j0 = Jstar[0]
    offs = sorted(j - j0 for j in Jstar)
    if offs[0] != 0:
        return False, None
    S = set(offs)
    k = len(S)
    if k & (k - 1) != 0:
        return False, None
    for a in offs:
        for b in offs:
            if (a ^ b) not in S:
                return False, None
    basis = []
    for x in sorted(offs):
        cur = x
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
    box = {0}
    for d in basis:
        box |= {x ^ d for x in box}
    return (box == S), sorted(basis)

def main():
    rows = []
    cnt_abs = Counter()
    odd_ge3 = []
    eights = []
    parity_violations = []
    box_failures = []
    for m in range(1, MMAX + 1):
        for lam in partitions(2*m):
            r = analyze(lam)
            a = r["absJstar"]
            cnt_abs[a] += 1
            if a >= 3 and a % 2 == 1:
                odd_ge3.append(r)
            if a == 8:
                eights.append(r)
            if a >= 2:
                if r["Jparity"] == "MIXED":
                    parity_violations.append(r)
                ok, gens = offsets_are_2adic_box(r["Jstar"])
                if not ok:
                    box_failures.append((r["lam"], r["Jstar"]))
            rows.append(r)
        print(f"  m={m:2d} done  ({sum(1 for _ in partitions(2*m))} shapes)", flush=True)

    print("\n=== |J*| census (all m<=%d) ===" % MMAX)
    for a in sorted(cnt_abs):
        print(f"  |J*|={a}: {cnt_abs[a]}")

    print("\n=== ODD |J*|>=3 (HEADLINE if nonempty) ===")
    if odd_ge3:
        print("  *** WALL BROKEN ***")
        for r in odd_ge3:
            print("   ", r["lam"], "J*=", r["Jstar"], "V=", r["V"])
    else:
        print("  NONE -- |J*| is even whenever >=2 (wall holds).")

    print("\n=== parity of J* members ===")
    if parity_violations:
        print("  MIXED parity found (would break a(s)/b(s) split):")
        for r in parity_violations:
            print("   ", r["lam"], r["Jstar"])
    else:
        print("  all tie shapes have single-parity J* (a(s)/b(s) split valid).")

    print("\n=== 2-adic box check on J* offsets ===")
    if box_failures:
        print(f"  {len(box_failures)} tie shapes whose offsets are NOT a 2-adic box:")
        for lam, J in box_failures[:20]:
            print("   ", lam, J)
    else:
        print("  every tie shape's J* offsets form an affine 2-adic box.")

    print("\n=== |J*|=8 cases ===")
    if eights:
        for r in eights:
            ok, gens = offsets_are_2adic_box(r["Jstar"])
            print("   ", r["lam"], "J*=", r["Jstar"], "box?", ok, "gens", gens)
    else:
        print("  none with m<=%d" % MMAX)

    # dump tie shapes for the Newton-polygon / involution scripts
    ties = [{"lam": r["lam"], "m": r["m"], "M": r["M"], "Jstar": r["Jstar"],
             "V": r["V"], "Jparity": r["Jparity"]} for r in rows if r["absJstar"] >= 2]
    with open("results-jstar-ties.json", "w") as f:
        json.dump(ties, f)
    print(f"\nwrote {len(ties)} tie shapes -> results-jstar-ties.json")

if __name__ == "__main__":
    main()
