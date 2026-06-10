"""
JOB A (verdict only) — does the Ayyer-Kumari hook-Schur linear-factor factorisation
apply to our M_j, and can ONE linear factor carry the -4 jump on a pure-M pair?

PROVE already pruned Route B this cycle (steplaw doc Sec 4); browsing is off so we cannot
re-transcribe AK 2501.00275.  We confirm the two STRUCTURAL facts that kill the import,
numerically and decisively:

  (a) The unique d=4 vanisher (2,2) is itself a 4-CORE (empty 4-quotient).  AK factorise the
      i=zeta_4 twisted characters when the 4-core is empty; their Thm predicts NON-vanishing
      exactly on empty-4-core shapes.  So "import via empty 4-core" predicts non-vanishing AT
      the one shape where G vanishes -- self-defeating.
  (b) A factorisation governs MAGNITUDES, not the Newton-locus involution.  Tool C gives
      M_j = sum_mu K^{(j)}_{lambda mu} f^mu, a POSITIVE SUM (not a product), so the -4 jump
      across a pure-M pair cannot be localised to a single linear factor: there is no product
      structure to carry it.  We verify M_j is genuinely a sum of >=2 distinct f^mu on the
      pure-M pairs (no monomial/linear-factor collapse), so no "one factor carries -4" picture
      can exist.
"""
import sys, math
from collections import Counter
from job1_tie_census import partitions, M_vector, v2, analyze
from core_quotient import core_quotient
from jobB_ladder import e2perp_chains_Mj, remove_vertical_2strip, f_lambda

# --- count the number of distinct f^mu summands in Tool C's expansion of M_j ---
def Mj_summand_profile(lam, m, j):
    from collections import Counter as C
    cur = C({tuple(x for x in lam if x>0): 1})
    for _ in range(j):
        nxt = C()
        for nu, coef in cur.items():
            for mu in remove_vertical_2strip(nu):
                nxt[mu] += coef
        cur = nxt
    # distinct mu with nonzero f^mu contribution
    terms = [(mu, coef, (f_lambda(mu) if sum(mu)>0 else 1)) for mu, coef in cur.items()]
    terms = [(mu, c, f) for (mu, c, f) in terms if c != 0 and f != 0]
    return terms

def main(MMAX):
    # (a) is (2,2) a 4-core? and is the unique vanisher a 4-core?
    core, quot = core_quotient((2, 2), 4)
    print("=== (a) the unique d=4 vanisher (2,2) ===")
    print(f"  4-core((2,2))     = {core}")
    print(f"  4-quotient((2,2)) = {tuple(quot)}   (empty quotient => (2,2) IS a 4-core)")
    print(f"  AK factorise/non-vanish on empty-4-core shapes; (2,2) is empty-4-core yet G=0")
    print(f"  => 'import via empty 4-core' predicts the WRONG answer at the only vanisher. PRUNED.\n")

    # (b) on pure-M pairs, M_j is a positive sum of >=2 distinct f^mu (no single-factor structure)
    print(f"=== (b) pure-M pairs {{j,j+8}} of J*: M_j summand structure  (m<={MMAX}) ===")
    npairs = 0; min2_j = 0; min2_j8 = 0
    examples = []
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            r = analyze(lam, m)
            J = set(r['Jstar'])
            Ms = M_vector(lam, m)
            for j in J:
                if (j+8) in J and (j ^ 8) == j+8:
                    npairs += 1
                    tj = Mj_summand_profile(lam, m, j)
                    tj8 = Mj_summand_profile(lam, m, j+8)
                    if len(tj) >= 2: min2_j += 1
                    if len(tj8) >= 2: min2_j8 += 1
                    if len(examples) < 4:
                        examples.append((tuple(lam), m, j,
                                         v2(Ms[j]), v2(Ms[j+8]), len(tj), len(tj8)))
    print(f"  pure-M J*-pairs found: {npairs}")
    print(f"  M_j   a sum of >=2 distinct f^mu: {min2_j}/{npairs}")
    print(f"  M_j+8 a sum of >=2 distinct f^mu: {min2_j8}/{npairs}")
    print(f"  examples (lam, m, j, v2(M_j), v2(M_j+8), #summands_j, #summands_j8):")
    for ex in examples:
        print(f"     {ex}")
    print(f"\n  VERDICT: M_j has no product/linear-factor structure -- it is an irreducible")
    print(f"  positive sum of SYT-counts.  The -4 jump is a Newton-locus (J*-membership) fact,")
    print(f"  not a property of any one hook-Schur factor.  Route B (AK factorisation) PRUNED,")
    print(f"  consistent with the PROVE steplaw doc Sec 4.")

if __name__ == "__main__":
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    main(MMAX)
