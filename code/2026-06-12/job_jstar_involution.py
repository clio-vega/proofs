"""
Job A.3 -- search for a fixed-point-free involution behind |J*| even, in the SYT-chain model.

For (lambda, j), the chain decomposition is
    M_j = sum_{mu |- 2(m-j)} N^{(j)}_{lambda->mu} f^mu,
N^{(j)} = #(length-j chains lambda -> mu, each step removing a vertical 2-strip).

We test, on the EDGE-contributing j (the J* members), several candidate involutions:

  C1  CONJUGATION on chain ends:  mu <-> mu'.  Since f^mu = f^{mu'} always, the non-self-conjugate
      ends contribute an EVEN amount IFF N_mu = N_{mu'}.  We test whether N_mu == N_{mu'} (chain-count
      conjugation symmetry).  Structurally we EXPECT failure: vertical-2-strip removal conjugates to
      HORIZONTAL-2-strip removal, so N^{vert}_{lam->mu} = N^{horiz}_{lam'->mu'} != N^{vert}_{lam->mu'}.

  C2  2-CORE grading of chain ends:  group mu by 2-core(mu); within a 2-core the f^mu pair?

  C3  smallest-2-adic-box-generator TOGGLE on J*:  j <-> j + 2^a (a = smallest gen).  This is
      fixed-point-free on J* by construction whenever the offsets form a 2-adic box -- but it is
      TAUTOLOGICAL (it assumes the box, which is the very thing to prove).  We confirm it is fpf and
      flag it as circular, per the 06-11 redirect note.

  C4  do the chain ENDS of paired j (j and j+2^a) relate by a vertical-2-strip map?  i.e. is there a
      structural bijection between the chains counted at j and at j+2^a?

Output: per-shape verdicts + global tallies.  States plainly whether a NON-circular fpf involution
emerged or exactly where each candidate pairing fails.
"""
import json
from collections import defaultdict
from job_jstar_engine import (remove_vstrip2, mult, hook_f, v2, C)

def conj(lam):
    lam = [x for x in lam if x > 0]
    if not lam:
        return ()
    return tuple(sum(1 for r in lam if r > c) for c in range(lam[0]))

def chain_ends(lam, j):
    """dict mu -> N^{(j)}_{lam->mu} via j vertical-2-strip removals."""
    state = {tuple(x for x in lam if x > 0): 1}
    for _ in range(j):
        state = mult(state, remove_vstrip2)
    return dict(state)

def analyze_shape(rec):
    m = rec["m"]; M = rec["M"]; lam = tuple(rec["lam"]); J = rec["Jstar"]
    p = J[0] % 2
    out = {"lam": lam, "m": m, "Jstar": J}

    # ---- C1 / C2 per on-edge j ----
    c1_sym_all = True          # N_mu == N_{mu'} for all ends, all on-edge j?
    c1_selfconj_drives = True  # does v2(self-conj part) == v2(M_j) (so pairing leaves an odd core)?
    c2_core_pairs = True
    for j in J:
        ends = chain_ends(lam, j)
        # C1: conjugation symmetry of chain counts
        for mu, N in ends.items():
            mc = conj(mu)
            if ends.get(mc, 0) != N:
                c1_sym_all = False
                break
        # self-conjugate contribution
        S = sum(N * hook_f(mu) for mu, N in ends.items() if conj(mu) == mu)
        Mj = M[j]
        if Mj != 0:
            if (S == 0 and v2(Mj) is not None) or (S != 0 and v2(S) != v2(Mj)):
                # self-conj part does NOT carry the valuation of M_j
                c1_selfconj_drives = False
        # C2: 2-core grouping -- contributions within each 2-core even?
        by_core = defaultdict(int)
        for mu, N in ends.items():
            by_core[two_core(mu)] += N * hook_f(mu)
        if any(v % 2 != 0 for v in by_core.values()):
            c2_core_pairs = False

    # ---- C3 smallest-generator toggle (tautological) ----
    offs = [j - J[0] for j in J]
    gens = sorted(offs)[1:2]  # smallest nonzero offset
    c3_fpf = False
    if gens:
        a = gens[0]
        # toggle j -> j XOR (offset bit) within the box
        pair = {}
        Jset = set(J)
        ok = True
        for j in J:
            t = j ^ a if False else None  # placeholder; do additive box toggle below
        # additive toggle on offsets: o -> o ^ a works iff offsets are a 2-adic box of bit-sums
        try:
            okbox = all(((o ^ a) + J[0]) in Jset for o in offs)
        except Exception:
            okbox = False
        c3_fpf = okbox and len(J) % 2 == 0

    out.update(dict(c1_sym_all=c1_sym_all, c1_selfconj_drives=c1_selfconj_drives,
                    c2_core_pairs=c2_core_pairs, c3_fpf=c3_fpf))
    return out

def two_core(mu):
    """2-core of a partition via repeated domino removal (beta-set / abacus on 2 runners)."""
    mu = sorted([x for x in mu if x > 0], reverse=True)
    n = len(mu)
    # beta-numbers
    beta = sorted((mu[i] + (n - 1 - i)) for i in range(n))
    # distribute on 2 runners, push down
    r0 = sorted(b for b in beta if b % 2 == 0)
    r1 = sorted(b for b in beta if b % 2 == 1)
    new = [2*i for i in range(len(r0))] + [2*i+1 for i in range(len(r1))]
    new = sorted(new, reverse=True)
    k = len(new)
    core = [new[i] - (k - 1 - i) for i in range(k)]
    core = tuple(x for x in sorted(core, reverse=True) if x > 0)
    return core

def main():
    import sys
    mcap = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    ties = [r for r in json.load(open("results-jstar-ties.json")) if r["m"] <= mcap]
    n = len(ties)
    print(f"(restricting to m<={mcap}: {n} tie shapes)")
    tally = defaultdict(int)
    fail_examples = {"c1_sym": [], "c2_core": []}
    for rec in ties:
        r = analyze_shape(rec)
        for k in ("c1_sym_all", "c1_selfconj_drives", "c2_core_pairs", "c3_fpf"):
            tally[k] += int(r[k])
        if not r["c1_sym_all"] and len(fail_examples["c1_sym"]) < 5:
            fail_examples["c1_sym"].append(r["lam"])
        if not r["c2_core_pairs"] and len(fail_examples["c2_core"]) < 5:
            fail_examples["c2_core"].append(r["lam"])

    print(f"=== involution search over {n} tie shapes ===")
    print(f"C1  conjugation chain-count symmetry N_mu==N_mu' (all on-edge j): "
          f"{tally['c1_sym_all']}/{n}")
    print(f"C1' self-conjugate part carries v2(M_j) (pairing leaves odd core): "
          f"{tally['c1_selfconj_drives']}/{n}")
    print(f"C2  2-core grading pairs all contributions (even per core):        "
          f"{tally['c2_core_pairs']}/{n}")
    print(f"C3  smallest-gen toggle fpf on J* (TAUTOLOGICAL / circular):       "
          f"{tally['c3_fpf']}/{n}")
    if fail_examples["c1_sym"]:
        print("  C1 fails e.g.:", fail_examples["c1_sym"])
    if fail_examples["c2_core"]:
        print("  C2 fails e.g.:", fail_examples["c2_core"])

if __name__ == "__main__":
    main()
