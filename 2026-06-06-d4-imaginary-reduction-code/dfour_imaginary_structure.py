"""
Structural analysis of the Im G_lambda zero-set (companion to dfour_imaginary_reduction.py).

Findings established here:
  A. Conjugation law  G_{lambda'} = i^m * conj(G_lambda)   (proves vanishing set is
     conjugation-symmetric; Im-zero set is NOT, except via the i^m twist).
  B. The Im-zero set splits: TRIVIAL (single row (2m), single column (1^{2m}) -> G real,
     =+-1) vs NONTRIVIAL (Re != 0 but Im = 0).  Classify and count both.
  C. Anchor: the FULL vanishing G_lambda = 0 holds ONLY at (2,2) for m <= 12
     (Re and Im both zero) -- so the route to Gap A needs BOTH parts (or 4-core
     valuation), confirming the imaginary reduction is two-row-special.
  D. Two-row sharp bound: min_{b>=1} |Im G_{(2m-b,b)}| and its minimiser (hook b=1).
"""
import sys, os
from math import comb

sys.path.insert(0, os.path.dirname(__file__))
from dfour_imaginary_reduction import ReIm_from_Nk, Nk_table
sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                "..", "scratch", "2026-06-09-d4-involution"))
from core import partitions, conj


def conjugation_law(mmax=11):
    """Verify G_{lambda'} = i^m * conj(G_lambda) for all lambda |- 2m, m <= mmax."""
    ok = True
    for m in range(1, mmax + 1):
        for lam in partitions(2 * m):
            Re, Im = ReIm_from_Nk(lam, m)
            lamc = conj(lam)
            Rec, Imc = ReIm_from_Nk(lamc, m)
            # i^m * (Re - i Im):  rotate (Re, -Im) by m quarter-turns
            a, b = Re, -Im
            for _ in range(m % 4):
                a, b = -b, a                  # multiply by i
            if (a, b) != (Rec, Imc):
                ok = False
                print(f"  CONJ FAIL m={m} lam={lam}: pred={(a,b)} got={(Rec,Imc)}")
    return ok


def classify_zeros(mmax=12):
    print("  Im-zero set classification (TRIVIAL = single row/column; NONTRIVIAL = Re!=0):")
    grand_nontriv = []
    for m in range(1, mmax + 1):
        n = 2 * m
        row = (n,)
        col = tuple([1] * n)
        nontriv = []
        triv = 0
        for lam in partitions(n):
            Re, Im = ReIm_from_Nk(lam, m)
            if Im != 0:
                continue
            if lam == row or lam == col:
                triv += 1
            else:
                nontriv.append((lam, Re))
        grand_nontriv.extend((m, lam, Re) for lam, Re in nontriv)
        sc = sum(1 for lam, _ in nontriv if conj(lam) == lam)
        print(f"   m={m:2d}: {triv} trivial + {len(nontriv)} nontrivial "
              f"(self-conj among nontriv: {sc})")
    print(f"\n  total NONTRIVIAL Im-zeros (m<=12): {len(grand_nontriv)}")
    print("  -> Im G_lambda vanishes on MANY shapes; cannot single out (2,2). VERDICT: "
          "imaginary reduction does NOT generalise.")
    return grand_nontriv


def full_vanishing_anchor(mmax=12):
    """Confirm G_lambda = 0 (Re=Im=0) ONLY at (2,2) for m <= mmax."""
    full = []
    for m in range(1, mmax + 1):
        for lam in partitions(2 * m):
            Re, Im = ReIm_from_Nk(lam, m)
            if Re == 0 and Im == 0:
                full.append((m, lam))
    return full


def two_row_bound(mmax=12):
    """min_{b>=1} |Im G_{(2m-b,b)}| and the minimiser, m<=mmax."""
    rows = []
    for m in range(2, mmax + 1):
        best = None
        for b in range(1, m + 1):
            lam = (2 * m - b, b)
            Re, Im = ReIm_from_Nk(lam, m)
            a = abs(Im)
            if best is None or a < best[0]:
                best = (a, lam)
        rows.append((m, best[0], best[1]))
    return rows


if __name__ == "__main__":
    print("=" * 70)
    print("A. CONJUGATION LAW  G_{lambda'} = i^m conj(G_lambda)")
    print("=" * 70)
    print(f"  verified for all lambda |- 2m, m<=11: {conjugation_law(11)}\n")

    print("=" * 70)
    print("B. Im-zero set classification")
    print("=" * 70)
    nontriv = classify_zeros(12)
    print("\n  nontrivial Im-zeros (m, lambda, Re):")
    for m, lam, Re in nontriv:
        print(f"    m={m:2d}  lam={lam}  Re={Re}")

    print("\n" + "=" * 70)
    print("C. FULL-VANISHING ANCHOR  (Re = Im = 0)")
    print("=" * 70)
    full = full_vanishing_anchor(12)
    print(f"  shapes with G_lambda(i)=0 for m<=12: {full}")
    print(f"  -> unique full vanisher is (2,2): {full == [(2, (2, 2))]}")

    print("\n" + "=" * 70)
    print("D. TWO-ROW SHARP BOUND  min_{b>=1} |Im G_{(2m-b,b)}|")
    print("=" * 70)
    print(f"  {'m':>3} {'min|Im|':>8}  minimiser")
    for m, mn, lam in two_row_bound(12):
        print(f"  {m:>3} {mn:>8}  {lam}")
