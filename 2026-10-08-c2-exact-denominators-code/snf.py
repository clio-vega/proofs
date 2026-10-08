"""
MAIN TARGET.  Claim:  P^(L) := Z-span{p_rho : l(rho)<=L}  =  Z-span{tilde m_nu : l(nu)<=L},
where tilde m_nu = d_nu m_nu.  Equivalently M^(L) = Mt^(L) . D^(L) with Mt unipotent over Z.
Consequence: coker(M^(L)) = (+)_{l(rho)<=L} Z/d_rho, i.e. SNF(M^(L)) = SNF(diag(d_rho)).
"""
from fractions import Fraction
from lattice import *
from sympy import Matrix, Integer
from sympy.matrices.normalforms import smith_normal_form
import sys

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8

def Mt_entry(lam, nu):
    """#{pi in Pi_{l(lam)} : lam^pi = nu}"""
    return sum(1 for pi in set_partitions(len(lam)) if coarsen(lam, pi) == nu)

def Mtinv_entry(lam, nu):
    """sum over the fibre of mu(hat0,pi)"""
    return sum(mu_hat0(pi) for pi in set_partitions(len(lam)) if coarsen(lam, pi) == nu)

print("=== (1) factorisation  M = Mt . D  entrywise, and unipotence of Mt ===")
fac_bad = []; uni_bad = []; nfac = 0; nfac_offdiag_nonzero = 0
for n in range(1, NMAX + 1):
    parts = partitions(n)
    for lam in parts:
        rowB = M_row_B(lam, n)
        for nu in parts:
            lhs = rowB.get(nu, 0)
            rhs = Mt_entry(lam, nu) * d_of(nu)
            nfac += 1
            if lhs != rhs: fac_bad.append((n, lam, nu, lhs, rhs))
            if lam != nu and Mt_entry(lam, nu) != 0: nfac_offdiag_nonzero += 1
        if Mt_entry(lam, lam) != 1: uni_bad.append(('diag', n, lam, Mt_entry(lam, lam)))
        # triangularity for length
        for nu in parts:
            if Mt_entry(lam, nu) != 0 and len(nu) > len(lam):
                uni_bad.append(('tri', n, lam, nu))
print("  entries checked            : %d   (n <= %d)" % (nfac, NMAX))
print("  NONZERO OFF-DIAGONAL in Mt : %d   <-- unipotence is non-vacuous only where these are" % nfac_offdiag_nonzero)
print("  factorisation failures     : %d" % len(fac_bad))
for b in fac_bad[:8]: print("   ", b)
print("  unipotence/triangularity failures: %d" % len(uni_bad))

print()
print("=== (2) Mt^{-1} is integral, and equals the Mobius fibre sums ===")
inv_bad = []; inv_checked = 0; inv_nonzero_offdiag = 0
for n in range(1, min(NMAX,7)+1):
    parts = partitions(n)
    Mt = Matrix([[Integer(Mt_entry(l, v)) for v in parts] for l in parts])
    Mtinv = Mt.inv()
    for i, l in enumerate(parts):
        for j, v in enumerate(parts):
            got = Mtinv[i, j]
            want = Mtinv_entry(l, v)
            inv_checked += 1
            if got != want: inv_bad.append((n, l, v, got, want))
            if got != 0:
                assert got == int(got), ("NONINTEGRAL", n, l, v, got)
                if i != j: inv_nonzero_offdiag += 1
print("  entries checked                 : %d" % inv_checked)
print("  nonzero off-diagonal in Mt^{-1} : %d" % inv_nonzero_offdiag)
print("  mismatches vs Mobius fibre sums : %d" % len(inv_bad))
print("  all entries integral            : yes (asserted entrywise)")

print()
print("=== (3) COKERNEL:  SNF(M^(L)) == SNF(diag(d_rho)) ? ===")
rows = []
snf_bad = []
noncyclic = 0; nontrivial = 0; total_nL = 0
for n in range(1, NMAX + 1):
    parts_all = partitions(n)
    for L in range(1, n + 1):
        parts = [r for r in parts_all if len(r) <= L]
        if not parts: continue
        M = Matrix([[Integer(M_row_B(l, n).get(v, 0)) for v in parts] for l in parts])
        D = Matrix.diag(*[Integer(d_of(r)) for r in parts])
        sM = smith_normal_form(M); sD = smith_normal_form(D)
        divM = [sM[i, i] for i in range(len(parts))]
        divD = [sD[i, i] for i in range(len(parts))]
        total_nL += 1
        if divM != divD: snf_bad.append((n, L, divM, divD))
        inv = [int(x) for x in divD if x != 1]
        if len(inv) >= 1: nontrivial += 1
        if len(inv) >= 2: noncyclic += 1
        if L == n:
            rows.append((n, [int(x) for x in divD if x != 1], int(D.det())))
print("  (n,L) pairs tested            : %d   (n <= %d)" % (total_nL, NMAX))
print("  pairs with NONTRIVIAL cokernel: %d   <-- the falsifiable ones (cokernel=0 is vacuous)" % nontrivial)
print("  pairs with NON-CYCLIC cokernel: %d   <-- where the group, not just the index, has content" % noncyclic)
print("  SNF mismatches                : %d" % len(snf_bad))
for b in snf_bad[:8]: print("   ", b)
print()
print("  L = n  (full matrix): invariant factors > 1 of Lambda_n^Z / P_n, and the index")
for n, inv, det in rows:
    print("    n=%d  index=%-8d invariant factors: %s" % (n, det, inv if inv else "(trivial)"))
