"""
JOB B  --  Map the SHARP (1+i)-adic Newton polygon of A_lambda(x) at x=i-1,
            i.e. the j-picture J* = argmin_j  val(j),  val(j) = j + 2 v2(C(m,j) M_j),
            for all lambda |- 2m, m <= 8.  Gather what a prove session needs for H1/H2.

Findings to establish/refine:
  H1 (parity class).  *** TRIVIAL ***:  val(j) = j + 2 v2(...) == j  (mod 2),
      so any two indices with equal val are congruent mod 2.  J* sits in one
      parity class automatically; "odd step" val(j+1)-val(j) is odd by construction.
      We verify the identity val(j)==j (mod 2) holds with no exception (sanity),
      then record H1 as a one-line corollary -- NOT a deep fact.

  H2 (2-adic box).  J* - j0 = { sum_{a in S} 2^a } for some generator set S of
      EVEN powers (forced by H1).  |J*| = 2^|S| in {1,2,4,8}.  Confirm to m=8;
      tabulate the generator multiset S and the |J*| distribution.

  S vs quotient.  Tabulate S against the 4-core and 4-quotient of lambda.
      Is S (the box generators) a function of the 4-quotient?  (06-11 said the
      survivor DEPTH is not; S may be coarser.)

  Anchor.  M_0 = f^lambda; tabulate v2(f^lambda) vs min_j val(j); when is 0 in J*?
"""
import math, csv
from collections import Counter, defaultdict

from job1_tie_census import partitions, v2, M_vector, analyze
from job_box_detail import is_affine_box
from quotient import t_quotient, f_syt


def val_vector(lam, m):
    Ms = M_vector(lam, m)
    val = {}
    for j in range(m + 1):
        c = math.comb(m, j) * Ms[j]
        if c != 0:
            val[j] = j + 2 * v2(c)
    return val, Ms


# ---------------------------------------------------------------- H1 triviality
def check_H1(MMAX=8):
    print("=" * 72)
    print("H1.  val(j) == j (mod 2) for every finite j  =>  J* single parity (TRIVIAL)")
    viol = 0
    tot = 0
    for m in range(1, MMAX + 1):
        for lam in partitions(2 * m):
            val, _ = val_vector(lam, m)
            for j, v in val.items():
                tot += 1
                if (v - j) % 2 != 0:
                    viol += 1
    print(f"   checked {tot} (shape,j) pairs (m<=8); val(j)!=j mod2 violations: {viol}")
    print("   => H1 holds with ZERO exceptions and is a one-line consequence of the")
    print("      definition val=j+2v2(.). No deep parity argument is needed.")
    return viol == 0


# ---------------------------------------------------------------- H2 box + S
def check_H2(MMAX=8):
    print("=" * 72)
    print("H2.  J* is an affine 2-adic box; tabulate generators S and |J*|.")
    box_ok = 0
    tie_tot = 0
    all_tot = 0
    size_dist = Counter()
    gen_dist = Counter()          # multiset of generators (offsets) over ties
    bad = []
    rows = []
    for m in range(1, MMAX + 1):
        for lam in partitions(2 * m):
            val, Ms = val_vector(lam, m)
            mu = min(val.values())
            J = sorted(k for k in val if val[k] == mu)
            all_tot += 1
            size_dist[len(J)] += 1
            if len(J) < 2:
                continue
            tie_tot += 1
            ok, gens = is_affine_box(J)
            if ok:
                box_ok += 1
                gen_dist[tuple(gens)] += 1
            else:
                bad.append((m, lam, J))
            core4, quot4 = t_quotient(lam, 4)
            rows.append(dict(m=m, lam="|".join(map(str, lam)), Jstar="|".join(map(str, J)),
                             Jsize=len(J), j0=J[0], gens="|".join(map(str, gens)) if ok else "NONBOX",
                             mu=mu, parity=J[0] % 2,
                             core4="|".join(map(str, core4)),
                             quot4=";".join("|".join(map(str, q)) for q in quot4)))
    print(f"   all shapes m<=8: {all_tot};  ties (|J*|>=2): {tie_tot}")
    print(f"   J* is an affine 2-adic box: {box_ok}/{tie_tot}")
    if bad:
        print("   NON-BOX ties:", bad[:10])
    print(f"   |J*| distribution (all shapes): {dict(sorted(size_dist.items()))}")
    print(f"   generator-set S distribution over ties (offset of each gen):")
    for g, c in sorted(gen_dist.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"      S(offsets)={g!s:14s} -> |J*|={2**len(g)} : {c} ties")
    return rows, box_ok, tie_tot


# ---------------------------------------------------------------- S vs 4-quotient
def S_vs_quotient(rows):
    print("=" * 72)
    print("S vs 4-quotient:  is the generator set a function of the 4-quotient?")
    # group ties by 4-quotient, see if S is constant within group
    by_quot = defaultdict(set)
    by_quot_core = defaultdict(set)
    for r in rows:
        if r['gens'] == "NONBOX":
            continue
        by_quot[r['quot4']].add(r['gens'])
        by_quot_core[(r['core4'], r['quot4'])].add(r['gens'])
    multi = sum(1 for v in by_quot.values() if len(v) > 1)
    print(f"   distinct 4-quotients among ties: {len(by_quot)}; "
          f"with >1 distinct S: {multi}")
    if multi:
        print("   => S is NOT a function of the 4-quotient alone. Examples:")
        shown = 0
        for q, ss in by_quot.items():
            if len(ss) > 1 and shown < 5:
                print(f"      quot4={q}  -> S in {sorted(ss)}")
                shown += 1
    # NB: (4-core,4-quotient) <-> lambda is a BIJECTION, so "S a function of
    # (core,quotient)" is vacuous (= S a function of lambda).  Do NOT report it.
    print("   NB: 4-core+4-quotient recover lambda bijectively, so that pair")
    print("       'determining' S is vacuous; the real test is the quotient alone (above).")


# ---------------------------------------------------------------- anchor
def anchor(MMAX=8):
    print("=" * 72)
    print("Anchor:  M_0 = f^lambda; v2(f^lambda) vs min val; when is 0 in J*?")
    f_eq = 0
    tot = 0
    zero_in_J = 0
    tie_zero_in_J = 0
    tie_tot = 0
    for m in range(1, MMAX + 1):
        for lam in partitions(2 * m):
            Ms = M_vector(lam, m)
            tot += 1
            if Ms[0] == f_syt(lam):
                f_eq += 1
            val, _ = val_vector(lam, m)
            mu = min(val.values())
            J = [k for k in val if val[k] == mu]
            if 0 in J:
                zero_in_J += 1
            if len(J) >= 2:
                tie_tot += 1
                if 0 in J:
                    tie_zero_in_J += 1
    print(f"   M_0 == f^lambda: {f_eq}/{tot}")
    print(f"   0 in J* (all shapes): {zero_in_J}/{tot}")
    print(f"   0 in J* (ties only):  {tie_zero_in_J}/{tie_tot}")


if __name__ == "__main__":
    check_H1(8)
    rows, box_ok, tie_tot = check_H2(8)
    S_vs_quotient(rows)
    anchor(8)
    with open("results/jobB-sharp-Jstar.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["m", "lam", "Jstar", "Jsize", "j0",
                                           "gens", "mu", "parity", "core4", "quot4"])
        w.writeheader()
        w.writerows(rows)
    print("=" * 72)
    print(f"wrote results/jobB-sharp-Jstar.csv ({len(rows)} tie shapes m<=8)")
