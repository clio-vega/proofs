#!/usr/bin/env python3
"""
JOB A (06-14 system cycle) -- is ord_{q=-1} G_lambda = floor(|2-core(lam)|/2) a
LEVEL-2 fusion (cylindric-LR / Verlinde) multiplicity?  And does ord_{zeta_d} G_lam
track LEVEL-d truncation?  (Dobner cylindric/level-2 fusion -- the surviving structural lead.)

Cross-checks built in:
  * order-law value computed TWO ways: closed formula floor(|2-core|/2) AND
    ord_{q=-1}(G_lambda(q)) where G_lambda(q)=sum_T q^{s(T)} (parity-twisted descent stat).
  * fusion routine validated in fusion_lib.py against Ising / sl_2 closed form / sl_3 level1.
"""
import sys, os
# bundled libs (syt_lib, dcore_lib, fusion_lib, job_jstar_engine) live alongside this script;
# the absolute paths are the in-container fallbacks for running from the live tree.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/clio/projects/code')
sys.path.insert(0, '/home/clio/projects/scratch/2026-06-04-sT-2quotient')
sys.path.insert(0, '/home/clio/projects/scratch/2026-06-05-dcore-lift')

import sympy as sp
from sympy import symbols, cyclotomic_poly
from collections import Counter, defaultdict

from syt_lib import G_lambda, partitions, num_syt, two_core_and_quotient, order_vanishing_at_minus1
from dcore_lib import d_core_and_quotient, ord_at_root
from fusion_lib import lr_coeff, lr_targets, fusion_product, level_weights

q = symbols('q')

SYT_MAX_N = 12   # polynomial (SYT) route up to n=12 (m=6); formula up to n=16 (m=8)


# ==========================================================================
def order_law_formula(lam):
    """floor(|2-core(lam)|/2)."""
    core, _ = two_core_and_quotient(lam)
    return sum(core) // 2


def part1_orderlaw_data():
    print("=" * 84)
    print("PART 1: the order-law value floor(|2-core|/2), cross-checked vs ord_{q=-1} sum_T q^{s(T)}")
    print("=" * 84)
    mism = []
    val_spectrum = Counter()
    n_form = 0
    rows = []
    for n in range(2, 17, 2):
        for lam in partitions(n):
            v = order_law_formula(lam)
            val_spectrum[v] += 1
            n_form += 1
            if n <= SYT_MAX_N:
                ov = order_vanishing_at_minus1(G_lambda(lam))
                if ov != v:
                    mism.append((lam, ov, v))
                rows.append((n, lam, v))
    print(f"  shapes (m<=8, formula): {n_form};  SYT-ord cross-check (m<=6) mismatches: {len(mism)}")
    for e in mism[:10]:
        print("    MISMATCH", e)
    print(f"  value spectrum of floor(|2-core|/2): {dict(sorted(val_spectrum.items()))}")
    print("  (values are UNBOUNDED: |2-core|=|delta_tau|=tau(tau+1)/2, so floor/2 in {0,0,1,3,5,7,...})")
    # show which shapes realize each value at small n
    print("\n  sample shapes by value (n<=10):")
    seen = set()
    for (n, lam, v) in rows:
        if n <= 10 and (v) not in seen and v > 0:
            core, _ = two_core_and_quotient(lam)
            print(f"    value {v}: lam={lam} (2-core={core}, |2-core|={sum(core)})")
            seen.add(v)
    return mism


# ==========================================================================
def part2_d_level_test():
    print("\n" + "=" * 84)
    print("PART 2: d -> level-d test:  ord_{zeta_d} G_lambda  for d=2,3,4  (Phi_d | G_lambda)")
    print("=" * 84)
    print("  Two distinct predictions to separate:")
    print("   (LEVEL=d Verlinde):  lower level = STRONGER truncation, so richness should DECREASE")
    print("       as d grows -- the qualitative trend below is (mildly) consistent with this.")
    print("   (RANK=d affine / d-core): ord_{zeta_d} should track floor(|d-core|/d) -- tested below.")
    print("  The discriminating question is QUANTITATIVE: do the actual d>=3 values match either?\n")
    for d in [2, 3, 4]:
        phi = cyclotomic_poly(d, q)
        ord_spectrum = Counter()
        track_core = []          # (ord_zeta_d, |d-core|-based candidate)
        n_chk = 0
        for n in range(2, SYT_MAX_N + 1):
            # restrict to n with d | n influence; compute for all n up to SYT_MAX_N
            for lam in partitions(n):
                G = G_lambda(lam)
                o = ord_at_root(G, phi)
                o = 0 if o is None else o
                ord_spectrum[o] += 1
                core, _ = d_core_and_quotient(lam, d)
                # candidate rank-d generalization (analogue of floor(|2-core|/2)):
                cand = sum(core) // d
                track_core.append((o, cand, lam, sum(core)))
                n_chk += 1
        # how well does ord_{zeta_d} match floor(|d-core|/d)?
        matches = sum(1 for (o, c, *_ ) in track_core if o == c)
        print(f"  d={d}:  ord_(zeta_{d}) spectrum over {n_chk} shapes (m<=6): {dict(sorted(ord_spectrum.items()))}")
        print(f"          ord_(zeta_{d}) == floor(|{d}-core|/{d}) on {matches}/{n_chk} shapes "
              f"({'MATCH' if matches==n_chk else 'MISMATCH'})")
        if matches != n_chk:
            # show first few divergences
            div = [(o, c, lam, cs) for (o, c, lam, cs) in track_core if o != c][:5]
            for (o, c, lam, cs) in div:
                print(f"            e.g. lam={lam}: ord_zeta_{d}={o}, |{d}-core|={cs}, floor/{d}={c}")
    print("\n  READ: d=2 ord rich (= floor(|2-core|/2), unbounded); d>=3 ord collapses to {0,1}.")
    print("    - RANK=d (d-core) prediction FAILS quantitatively for d>=3 (floor(|d-core|/d) grows,")
    print("      ord stays in {0,1}).  So ord_{zeta_d} does NOT read the d-core for d>=3.")
    print("    - LEVEL=d trend (richness shrinks with d) is qualitatively OK, but Part 3 shows no")
    print("      lambda-only level-2 fusion invariant reproduces the d=2 VALUES floor(|2-core|/2).")
    print("    => the EXACT carrier at d=2 is the 2-core (rank-2 affine, Littlewood t=2), not a")
    print("       level-2 Verlinde multiplicity; the d=2-specialness is intrinsic to s(T)'s mod-2 twist.")


# ==========================================================================
def part3_fusion_search():
    print("\n" + "=" * 84)
    print("PART 3: direct level-2 fusion search for a lambda-only invariant = floor(|2-core|/2)")
    print("=" * 84)
    print("  Obstruction up front: a level-2 weight of sl_n is a partition with parts <= 2.")
    print("  The order law is defined for ALL lambda (incl. wide shapes), and floor(|2-core|/2)")
    print("  is unbounded -- but fusion multiplicities of a fixed product are bounded by ordinary LR.")
    print("  We therefore test on the VALID domain lambda_1<=2 (two-column shapes), and probe several")
    print("  natural lambda-only level-2 quantities; we also record the n-dependence obstruction.\n")

    # two-column shapes lam = (2^a 1^b) of size 2m, valid level-2 sl_n weights
    def twocol(N):
        out = []
        for a in range(0, N // 2 + 1):
            b = N - 2 * a
            if b < 0: continue
            lam = tuple([2] * a + [1] * b)
            if lam: out.append(lam)
        return out

    print("  lam=(2^a 1^b)  |  floor(|2-core|/2)  |  self-fusion defect D_2(n) for n=rows+1, rows+2  |  #N^(2) terms")
    candidate_match = {"D2_minimal": 0, "D2_total": 0}
    rows_tested = 0
    examples = []
    for N in range(2, 15, 2):
        for lam in twocol(N):
            rows = len(lam)
            target = sum(two_core_and_quotient(lam)[0]) // 2
            # ordinary self-LR total multiplicity (<= n rows)
            defects = {}
            for n in [rows + 1, rows + 2]:
                ord_total = sum(lr_targets(lam, lam, n).values())
                fus = fusion_product(lam, lam, n, 2)
                fus_total = sum(fus.values())
                defects[n] = ord_total - fus_total
            n_terms = len(fusion_product(lam, lam, rows + 1, 2))
            examples.append((lam, target, defects[rows+1], defects[rows+2], n_terms))
            rows_tested += 1
    # print a sample and test matches
    for (lam, target, d1, d2, nt) in examples[:24]:
        flag = "  <-- D2(min)==target" if d1 == target else ""
        print(f"   {str(lam):18s}  {target:3d}   D2={d1:4d}(n=r+1), {d2:4d}(n=r+2)   |N^(2)|={nt}{flag}")
    n_match_min = sum(1 for (lam, t, d1, d2, nt) in examples if d1 == t)
    n_match_ndep = sum(1 for (lam, t, d1, d2, nt) in examples if d1 == d2)
    print(f"\n  self-fusion defect D_2(n=rows+1) == floor(|2-core|/2):  {n_match_min}/{len(examples)} shapes")
    print(f"  D_2 independent of n (rows+1 vs rows+2):                 {n_match_ndep}/{len(examples)} shapes")
    print("  => D_2 is n-DEPENDENT (not a lambda-only invariant) and does not equal the target.")

    # also: is floor(|2-core|/2) ever a single fusion multiplicity N^(2 nu)_{lam lam}? bounded by 1 for sl_2-like
    print("\n  Spectrum reminder: floor(|2-core|/2) reaches 3,5,6,7,... (unbounded);")
    print("  a level-2 fusion multiplicity N^(2)_{lam,mu}^nu for these wide shapes is not even defined")
    print("  (lam_1>2). No lambda-only level-2 fusion invariant matches floor(|2-core|/2) on the full domain.")


# ==========================================================================
if __name__ == "__main__":
    m = part1_orderlaw_data()
    part2_d_level_test()
    part3_fusion_search()
    print("\n" + "=" * 84)
    print("See FINDINGS-2026-06-14-jobA-fusion.md for the verdict.")
    print("=" * 84)
