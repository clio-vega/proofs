"""
TASK 2 (cont.) -- pin down the core x runner interaction Delta(core, runner).

task2_runner.py established:
  * runner-MULTISET does NOT resolve the 12 exceptions (K2 still 12 non-constant);
  * the per-runner BEAD-COUNT vector t_j of the full lambda abacus DOES (K3 -> 0).
  * at FIXED core, residual of a single-box quotient varies with the runner the box
    sits on, in a structured way:
        core=(5,2,1,1,1):   {0:4, 1:2, 2:2, 3:4}
        core=(6,3,1,1,1):   {0:2, 1:4, 2:2, 3:4}
   i.e. residual = base + Delta where Delta is high on two runners, low on two,
   and WHICH two depends on the core.

This script computes, for each 4-core c (canonical r = 4*ceil(ell/4)):
   - canonical bead-count vector  t_j(c)  and the bead positions on each runner,
   - the residual of a single box placed on runner j  (from m<=MMAX data),
and tests the conjecture that the "high-residual" runners are exactly those whose
canonical bead-count t_j(c) has a particular parity / value, i.e.
   Delta(c, j) = phi( t_j(c) )   for an explicit phi.
"""
import sys, os
sys.path.insert(0, "/home/clio/projects/scratch/2026-06-09-d4-involution")
from collections import defaultdict
from task2_runner import (build_rows, abacus_data_4, beta_set, hook_f, v2, vpi)

def canonical_core_beads(core):
    """Canonical 4-runner bead data of the core itself (r = 4*ceil(ell/4))."""
    ell = len([x for x in core if x > 0])
    r = ell
    while r % 4 != 0:
        r += 1
    if r == 0:
        r = 4
    betas = beta_set(core, r)
    runners = defaultdict(list)
    for b in betas:
        runners[b % 4].append(b // 4)
    t = tuple(len(runners[j]) for j in range(4))
    return r, t, {j: sorted(runners[j]) for j in range(4)}

def main():
    MMAX = 11
    rows = build_rows(MMAX)

    # single-box residual by (core, runner)
    sb = defaultdict(dict)   # core -> {runner: residual}
    for r in rows:
        if r['res'] is None or r['qsz'] != 1:
            continue
        runner = r['rms'][0]
        sb[r['core']].setdefault(runner, set()).add(r['res'])

    print("=" * 78)
    print("Single-box residual by runner, with canonical core bead-counts t_j(c)")
    print("=" * 78)
    records = []
    for core in sorted(sb, key=lambda c: (sum(c), c)):
        d = sb[core]
        if not all(len(v) == 1 for v in d.values()):
            continue
        resd = {j: next(iter(v)) for j, v in d.items()}
        rr, t, beads = canonical_core_beads(core)
        records.append((core, rr, t, resd))

    # Look for the rule:  is residual(core, runner=j) determined by (something canonical)?
    # Hypothesis A: residual depends on t_j(c) (count on that runner) + a base.
    print("\nHYP A: residual(c,j) as a function of t_j(c):")
    tableA = defaultdict(set)
    for core, rr, t, resd in records:
        for j, res in resd.items():
            tableA[t[j]].add(res)
    for tj in sorted(tableA):
        print(f"   t_j={tj}:  residuals={sorted(tableA[tj])}")

    # Hypothesis B: residual(c,j) - base(c) depends on t_j(c) - min_k t_k, etc.
    print("\nHYP B: residual(c,j) vs (t_j(c), parity of t_j, core size):")
    for core, rr, t, resd in records[:20]:
        line = "  ".join(f"r{j}:res={resd.get(j,'-')},t={t[j]}" for j in range(4))
        print(f"   core={str(core):22s} r={rr} t={t}  | {line}")

    # Hypothesis C: the HIGH-residual runners vs the bead-count vector
    print("\nHYP C: which runners are high-residual vs t_j ranking:")
    for core, rr, t, resd in records:
        if len(set(resd.values())) < 2:
            continue
        hi = max(resd.values()); lo = min(resd.values())
        hi_runners = sorted(j for j, v in resd.items() if v == hi)
        # relate to t_j: are high runners those with EXTREME t_j (max or min)?
        tmax = max(t); tmin = min(t)
        hi_by_t = sorted(j for j in range(4) if t[j] in (tmax, tmin))
        print(f"   core={str(core):22s} t={t}  hi_runners={hi_runners} (res {lo}->{hi})"
              f"  | runners with extreme t_j={hi_by_t}  match={hi_runners==hi_by_t}")

    # Hypothesis D: residual = base + 2*[number of core beads strictly below the
    # inserted bead on its runner is ODD]?  Try: high iff t_j(c) is EVEN / ODD.
    print("\nHYP D: high-residual runner  <=>  parity of t_j(c):")
    consistent_even = consistent_odd = True
    for core, rr, t, resd in records:
        if len(set(resd.values())) < 2:
            continue
        hi = max(resd.values())
        for j, v in resd.items():
            is_hi = (v == hi)
            if is_hi != (t[j] % 2 == 0):
                consistent_even = False
            if is_hi != (t[j] % 2 == 1):
                consistent_odd = False
    print(f"   high <=> t_j even : {consistent_even}")
    print(f"   high <=> t_j odd  : {consistent_odd}")

if __name__ == "__main__":
    main()
