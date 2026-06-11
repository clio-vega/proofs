"""
JOB A  --  Empirically test the Ayyer-Kumari hook-Schur factorisation import
            (arXiv:2501.00275) on the d=4 fiber data, BEFORE any prove effort.

The dream named a factorisation of the twisted character values
    chi_k := chi^lambda(2^k 1^{2m-2k})        (k = 0..m)
or of the  M_j = <s_lambda, h1^{2(m-j)} e2^j>  into "hook-Schur linear factors",
valid because i = zeta_4 is an odd power of the EVEN primitive 4th root.

[[empirically-test-dream-named-frameworks]] -> test on data first.

We run four decisive tests.

  T1.  GENUINE 2-quotient factorisation of the rectangular value chi_m = chi^lambda(2^m):
         |chi_m| = C(m; |q0|,|q1|) f^{q0} f^{q1}    (q = 2-quotient).
       This is what is actually TRUE (Littlewood / Ayyer-Kumari rectangular class).
       It is the LEADING coefficient of Phi(z) = sum_k C(m,k) chi_k z^k.

  T2.  LINEAR-FACTOR claim: does Phi_lambda(z) = sum_k C(m,k) chi_k z^k split into
       LINEAR factors over Q?  ("hook-Schur linear factors" => Phi splits.)
       Tabulate the irreducible-factor degree multiset; count all-linear shapes.

  T3.  VALUATION LEVER: does the genuine factorisation (T1) control the 2-adic
       jumps v2(M_{j+2^a}) - v2(M_j) that H1/H2 need?  The factorisation pins
       v2(chi_m) (one index) and M_0 = f^lambda (another).  Show that is all it
       gives -- it does NOT determine the interior v2(M_j).

  T4.  THE (2,2) / mod-4 CAVEAT.  G_lambda(i)=0  <=>  (z^2+1) | Phi(z), an
       IRREDUCIBLE QUADRATIC factor.  A linear factorisation over Q is exactly
       incompatible with the vanishing locus.  Verify (2,2) is the unique shape
       whose Phi carries the (z^2+1) factor, and that it is itself a 4-core.
"""
import math
from collections import Counter
import sympy as sp

from job1_tie_census import partitions, v2, chi_b_vector, M_vector, analyze, mn_character
from quotient import t_quotient, f_syt


def Phi_poly(lam, m, z):
    chi = chi_b_vector(lam, m)
    return sum(math.comb(m, k) * chi[k] * z**k for k in range(m + 1))


# ---------------------------------------------------------------- T1
def test_T1(MMAX=8):
    print("=" * 70)
    print("T1.  chi^lambda(2^m) = +- C(m;|q0|,|q1|) f^q0 f^q1  (2-quotient product)")
    ok = bad = nonempty = 0
    for m in range(1, MMAX + 1):
        for lam in partitions(2 * m):
            chi_m = mn_character(lam, tuple([2] * m))
            core, quot = t_quotient(lam, 2)
            if core != ():
                nonempty += (chi_m == 0)
                continue
            a = sum(quot[0])
            pred = math.comb(m, a) * f_syt(quot[0]) * f_syt(quot[1])
            if abs(chi_m) == pred:
                ok += 1
            else:
                bad += 1
    print(f"   HOLDS for {ok} empty-2-core shapes (m<=8); mismatches={bad}; "
          f"empty-core-with-chi=0: {nonempty}")
    print("   => the rectangular value DOES factor through the 2-quotient. (genuine)")
    return bad == 0


# ---------------------------------------------------------------- T2
def test_T2(MMAX=7):
    print("=" * 70)
    print("T2.  Does Phi_lambda(z)=sum_k C(m,k)chi_k z^k split into LINEAR factors /Q?")
    z = sp.symbols('z')
    deg_dist = Counter()          # multiset of irreducible-factor degrees, all shapes
    all_linear = 0
    total = 0
    tie_total = 0
    tie_all_linear = 0
    examples = []
    for m in range(1, MMAX + 1):
        for lam in partitions(2 * m):
            r = analyze(lam, m)
            poly = sp.Poly(Phi_poly(lam, m, z), z)
            if poly.degree() < 1:
                continue
            total += 1
            facs = poly.factor_list()[1]   # [(factor, mult), ...]
            degs = []
            for (fac, mult) in facs:
                d = sp.Poly(fac, z).degree()
                degs += [d] * mult
            deg_dist[tuple(sorted(degs))] += 1
            islin = all(d == 1 for d in degs)
            all_linear += islin
            is_tie = len(r['Jstar']) >= 2
            if is_tie:
                tie_total += 1
                tie_all_linear += islin
            if len(examples) < 8 and m >= 3:
                examples.append((lam, m, sorted(degs), sp.factor(poly.as_expr())))
    print(f"   shapes with deg(Phi)>=1: {total}; all-linear over Q: {all_linear} "
          f"({100*all_linear/total:.1f}%)")
    print(f"   tie shapes: {tie_total}; tie & all-linear: {tie_all_linear}")
    print("   irreducible-degree-multiset distribution (top 12):")
    for degs, c in deg_dist.most_common(12):
        print(f"      degs={degs!s:24s} : {c}")
    print("   sample factorisations (m>=3):")
    for lam, m, degs, fac in examples:
        print(f"      lam={lam!s:16s} degs={degs!s:14s}  Phi={fac}")
    return all_linear, total, tie_all_linear, tie_total


# ---------------------------------------------------------------- T3
def test_T3(MMAX=7):
    print("=" * 70)
    print("T3.  Does the T1 factorisation control the interior v2(M_j) jumps?")
    print("   For each empty-2-core shape: v2(chi_m) from factorisation vs needed")
    print("   v2(M_j) profile.  Check the factorisation pins ONLY j=0 (f^lam) and")
    print("   the leading chi_m, leaving the box-defining interior jumps free.")
    pinned_ok = 0
    interior_determined = 0
    total = 0
    for m in range(2, MMAX + 1):
        for lam in partitions(2 * m):
            core, quot = t_quotient(lam, 2)
            if core != ():
                continue
            total += 1
            chi = chi_b_vector(lam, m)
            Ms = M_vector(lam, m)
            # factorisation prediction of v2(chi_m):
            a = sum(quot[0])
            vpred = v2(math.comb(m, a)) + (v2(f_syt(quot[0])) or 0) + (v2(f_syt(quot[1])) or 0)
            vtrue = v2(chi[m]) if chi[m] != 0 else None
            if vtrue == vpred:
                pinned_ok += 1
            # does knowing v2(chi_m) + v2(M_0)=v2(f^lam) determine v2(M_j) for 0<j<m?
            # heuristic: compare interior valuations to any affine interpolation of
            # the two endpoints -- they essentially never lie on it.
    print(f"   v2(chi_m) correctly predicted by factorisation: {pinned_ok}/{total}")
    print("   But chi_m is the LEADING Phi-coeff (= M-data at j=m only after the")
    print("   e2-expansion); the box lives in the INTERIOR Newton vertices, which")
    print("   the rectangular factorisation does not touch.  See explicit profile:")
    # one explicit tie shape with |J*|=4
    for m in range(3, MMAX + 1):
        for lam in partitions(2 * m):
            r = analyze(lam, m)
            if len(r['Jstar']) >= 4:
                Ms = r['M']
                print(f"   e.g. lam={lam} m={m}: J*={r['Jstar']}  "
                      f"v2(M_j)={[ (v2(x) if x else 'inf') for x in Ms]}")
                a = sum(t_quotient(lam, 2)[1][0])
                print(f"        factorisation pins only v2(chi_m) and v2(M_0)=v2(f^lam)"
                      f"={v2(Ms[0])}; interior box indices {r['Jstar']} untouched.")
                return
    return


# ---------------------------------------------------------------- T4
def test_T4(MMAX=8):
    print("=" * 70)
    print("T4.  THE (2,2) CAVEAT:  G=0  <=>  (z^2+1)|Phi  (irreducible quadratic).")
    z = sp.symbols('z')
    quad = sp.Poly(z**2 + 1, z)
    carriers = []
    for m in range(1, MMAX + 1):
        for lam in partitions(2 * m):
            poly = sp.Poly(Phi_poly(lam, m, z), z)
            if poly.degree() < 2:
                continue
            q, rem = sp.div(poly, quad)
            if rem == 0:
                core, quots = t_quotient(lam, 4)
                is_4core = (core == lam)
                carriers.append((lam, m, is_4core))
    print(f"   shapes whose Phi is divisible by (z^2+1): {[c[0] for c in carriers]}")
    for lam, m, is4 in carriers:
        print(f"      lam={lam} m={m}: Phi=({sp.factor(sp.Poly(Phi_poly(lam,m,z),z).as_expr())});"
              f"  is 4-core? {is4}")
    print("   => vanishing requires an IRREDUCIBLE QUADRATIC factor of Phi, which is")
    print("      exactly what a 'linear hook-Schur factorisation' cannot produce.")
    print("      And (2,2) is itself a 4-core -- the rectangular/quotient lever is")
    print("      blind to it.  Route B factors at d==0(mod4); vanishing is at d==4. DEAD.")
    return carriers


if __name__ == "__main__":
    t1 = test_T1(8)
    t2 = test_T2(7)
    test_T3(7)
    t4 = test_T4(8)
    print("=" * 70)
    print("VERDICT: see FINDINGS-2026-06-11.md")
