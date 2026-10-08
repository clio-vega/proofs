"""
Q395 GATE, and the non-vacuity census that must precede any verdict.

Gate claim: M^(L)_{rho nu} = <p_rho,h_nu> equals the (truncated) m-to-p transition
matrix entry [m_nu] p_rho.  Engine C computes the left side with no reference to m;
engine B computes the right side with no reference to h or z.
"""
from fractions import Fraction
from lattice import *
import sys

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8

tot = 0; nonzero = 0; mismatch = []
falsifiable = 0
for n in range(1, NMAX + 1):
    parts = partitions(n)
    colC = {nu: M_col_C(nu) for nu in parts}
    for rho in parts:
        rowB = M_row_B(rho, n)
        for nu in parts:
            b = Fraction(rowB.get(nu, 0))
            c = colC[nu].get(rho, Fraction(0))
            tot += 1
            if b != 0 or c != 0: nonzero += 1
            # falsifiability: an entry is NOT forced by triangularity alone if it is
            # in the support, i.e. nu is a coarsening of rho.  Entries outside the
            # support are forced to 0 on both sides by the support theorem, so an
            # agreement there is rank-one.
            is_coarsening = any(coarsen(rho, pi) == nu for pi in set_partitions(len(rho)))
            if is_coarsening: falsifiable += 1
            if b != c: mismatch.append((n, rho, nu, b, c))

print("GATE  engine B ([m_nu]p_rho, monomial expansion) vs engine C (z_rho [p_rho]h_nu, Hall in p-basis)")
print("  n <= %d" % NMAX)
print("  entries compared        : %d" % tot)
print("  entries in the support  : %d   <-- the falsifiable ones" % falsifiable)
print("  entries outside support : %d   (forced 0 on both sides; agreement here is rank one)" % (tot - falsifiable))
print("  nonzero on either side  : %d" % nonzero)
print("  MISMATCHES              : %d" % len(mismatch))
for m in mismatch[:10]: print("   ", m)

# engine A as a third reading, on the falsifiable set only
mmA = 0; checkedA = 0
for n in range(1, min(NMAX, 7) + 1):
    for rho in partitions(n):
        rowB = M_row_B(rho, n)
        for nu in partitions(n):
            a = M_entry_A(rho, nu)
            b = rowB.get(nu, 0)
            checkedA += 1
            if a != b: mmA += 1
print("  engine A (function count) vs B: %d compared, %d mismatches (n<=%d)" % (checkedA, mmA, min(NMAX,7)))
