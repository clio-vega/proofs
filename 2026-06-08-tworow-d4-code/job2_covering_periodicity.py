"""
job2_covering_periodicity.py  --  Job 2.3 deepened: covering prime sets and b-periodicity.

For b == 2,3 (mod 4), q_b = primitive(Q_b) is IRREDUCIBLE (degree floor(b/2)), so its
roots mod p are exactly the roots of the primitive integer form of Q_b = I_b/(m)_{R+1}.
We therefore test "q_b root-free mod p" WITHOUT factoring (fast), over a long range of b,
to (a) test whether the covering prime sets persist and (b) detect b-periodicity of the
root-free property -- the only route to a UNIFORM (all-b) proof for these classes.

If, for a fixed prime p, {b in class : q_b root-free mod p} is a union of arithmetic
progressions in b, then a finite set of primes whose progressions cover the whole class
gives a uniform proof.  We hunt for exactly that.
"""

import sympy as sp
from sympy import Poly, primerange, div
from dfour_tworow import Ib, m

def Qb_primitive_modp_rootfree(b, p):
    """True iff primitive(Q_b) has no root in F_p.  No factoring."""
    P = Ib(b)
    R = (b - 1) // 2
    for r in range(R + 1):
        P = div(P, Poly(m - r, m, domain='QQ'), m)[0]
    # P is Q_b over QQ; clear denominators -> primitive integer poly
    Pz = Poly(P.as_expr(), m, domain='QQ')
    lcm_den = 1
    for c in Pz.all_coeffs():
        lcm_den = sp.ilcm(lcm_den, sp.Rational(c).q)
    Pint = Poly((P.as_expr() * lcm_den), m, domain='ZZ')
    cont = Pint.content()
    Pint = Poly(Pint.as_expr() / cont, m, domain='ZZ')
    fp = Poly(Pint.as_expr(), m, modulus=p)
    if fp.degree() < 1:
        return None  # degenerate mod p (p | leading) -- treat separately
    return not any(fp.eval(a) % p == 0 for a in range(p))

def rootfree_table(cls, bmax, primes):
    """Return dict p -> set of b in class (<=bmax) that are root-free mod p."""
    bs = [b for b in range(6, bmax + 1) if b % 4 == cls]
    table = {p: set() for p in primes}
    for b in bs:
        for p in primes:
            rf = Qb_primitive_modp_rootfree(b, p)
            if rf:
                table[p].add(b)
    return bs, table

def detect_period(bs_class, rootfree_b, cls, maxperiod=24):
    """Look for a period P (multiple of 4) s.t. root-free status depends only on b mod P."""
    for P in range(4, maxperiod + 1, 4):
        ok = True
        groups = {}
        for b in bs_class:
            groups.setdefault(b % P, []).append(b in rootfree_b)
        for res, vals in groups.items():
            if len(set(vals)) != 1:
                ok = False
                break
        if ok:
            rf_res = sorted(res for res, vals in groups.items() if vals[0])
            return P, rf_res
    return None, None

def main():
    print("=" * 78)
    print("Covering prime sets & b-periodicity for b == 2,3 (mod 4)")
    print("=" * 78)
    candidates = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
    BMAX = 62  # keep runtime bounded (deg up to ~30)
    for cls in (2, 3):
        print(f"\n--- class b == {cls} (mod 4), b in [6,{BMAX}] ---")
        bs, table = rootfree_table(cls, BMAX, candidates)
        print(f"  members: {bs}")
        # per-prime root-free sets + periodicity
        print(f"  {'p':>3} | {'#root-free':>10} | root-free b           | period (b mod P -> residues)")
        useful = []
        for p in candidates:
            rf = sorted(table[p])
            if not rf:
                continue
            P, res = detect_period(bs, set(rf), cls)
            perstr = f"P={P}: residues {res}" if P else "no period<=24"
            print(f"  {p:>3} | {len(rf):>10} | {str(rf):20s} | {perstr}")
            useful.append(p)
        # smallest covering set among candidates (each member root-free mod >=1)
        from itertools import combinations
        found = None
        for size in range(1, 5):
            for combo in combinations(useful, size):
                if all(any(b in table[p] for p in combo) for b in bs):
                    found = combo
                    break
            if found:
                break
        print(f"  => smallest covering set (tested members): {found} (size {len(found) if found else '>4'})")
        # if a covering set has periodic pieces, report the uniform decomposition
        if found:
            print("     periodic decomposition of the cover:")
            covered = set()
            for p in found:
                P, res = detect_period(bs, table[p], cls)
                tag = f"b mod {P} in {res}" if P else f"b in {sorted(table[p])}"
                print(f"       mod {p}: root-free for {tag}")

    print("\nNOTE: periodicity is over TESTED members only; a proven period would need a")
    print("recurrence argument (I_b mod p is governed by the order-? recurrence mod p).")


if __name__ == "__main__":
    main()
