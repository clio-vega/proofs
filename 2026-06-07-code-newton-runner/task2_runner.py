"""
TASK 2 -- the Gap-A abacus-runner interaction term.

Background (FINDINGS-4core-residual.md):
  residual(lambda) := v_pi(G_lambda) - v2(f^lambda) - v_pi(G_4core) + v2(f^4core).
  The near-formula  residual = f(v_pi(G_core), v2(f_core), 4-quotient)  held 775/787;
  at fixed core (5,2,1,1,1) a single quotient box gives residuals 2/4/4 depending on
  WHICH abacus runner it sits on.  Hypothesis: the runner position is the missing
  variable.

This script, for all lambda |- 2m, m <= MMAX:
  1. computes James-Kerber data on a 4-runner abacus with r = 4*ceil(ell/4) beads:
       - 4-core, ORDERED 4-quotient,
       - the runner index of each quotient box (= component index j),
       - runner-multiset (j with multiplicity |quotient^(j)|),
       - per-runner bead counts t_j of the FULL lambda abacus (the core's "charge"),
  2. computes residual,
  3. tests a LADDER of grouping keys (coarse -> fine) and reports #non-constant groups:
       K0: (vpc, v2fc, UNORDERED quotient shape)        -- core only via valuations
       K1: K0 + which-runner  == (vpc, v2fc, ORDERED quotient)
       K2: (vpc, v2fc, ORDERED quotient, runner-multiset)   (TASK 2 literal key)
       K3: (vpc, v2fc, ORDERED quotient, bead-count multiset {t_j})
       K4: (core, ORDERED quotient)                      -- complete invariant (sanity)
  4. zooms into the previous exception class and the (5,2,1,1,1) example,
  5. if a key makes residual a function, conjectures the correction term.

Cross-checks: core+4*|quot| = n; (core,quotient) is injective; unique vanisher (2,2).
"""
import sys, os
sys.path.insert(0, "/home/clio/projects/scratch/2026-06-09-d4-involution")
sys.path.insert(0, "/home/clio/projects/scratch/2026-06-07-prove")
from core import psi_power_schur
from math import factorial
from functools import lru_cache
from collections import defaultdict
import csv

RESULTS = "/home/clio/projects/code/results"
os.makedirs(RESULTS, exist_ok=True)

def v2(n):
    if n == 0:
        return float('inf')
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

def vpi(re, im):
    if re == 0 and im == 0:
        return float('inf')
    return v2(re * re + im * im)

def hook_f(lam):
    lam = [x for x in lam if x > 0]
    if not lam:
        return 1
    n = sum(lam)
    conj = [sum(1 for r in lam if r > c) for c in range(lam[0])]
    prod = 1
    for i, row in enumerate(lam):
        for j in range(row):
            arm = row - (j + 1)
            leg = conj[j] - (i + 1)
            prod *= (arm + leg + 1)
    return factorial(n) // prod

def beta_set(lam, r):
    lam = list(lam) + [0] * (r - len(lam))
    return [lam[i] + (r - 1 - i) for i in range(r)]

def beta_to_part(betas):
    betas = sorted(betas, reverse=True)
    r = len(betas)
    part = [betas[i] - (r - 1 - i) for i in range(r)]
    return tuple(x for x in part if x > 0)

def abacus_data_4(lam):
    """Return (core, ordered_quotient, runner_bead_counts t_j, runner_multiset)."""
    lam = [x for x in lam if x > 0]
    ell = len(lam)
    r = ell
    while r % 4 != 0:
        r += 1
    if r == 0:
        r = 4
    betas = beta_set(lam, r)
    runners = {j: [] for j in range(4)}
    for b in betas:
        runners[b % 4].append(b // 4)
    quotient = tuple(beta_to_part(runners[j]) for j in range(4))
    t = tuple(len(runners[j]) for j in range(4))           # bead counts per runner
    core_betas = []
    for j in range(4):
        for k in range(t[j]):
            core_betas.append(4 * k + j)
    core = beta_to_part(core_betas)
    # runner multiset: runner j with multiplicity |quotient^(j)|
    runner_ms = tuple(sorted(j for j in range(4) for _ in range(sum(quotient[j]))))
    return core, quotient, t, runner_ms

@lru_cache(maxsize=None)
def Gtab(m):
    return psi_power_schur(m)

def G_of(part):
    n = sum(part)
    if n == 0:
        return (1, 0)
    return Gtab(n // 2).get(tuple(part), (0, 0))

def build_rows(MMAX):
    rows = []
    for m in range(1, MMAX + 1):
        for lam, (re, im) in Gtab(m).items():
            vp = vpi(re, im)
            fl = hook_f(lam); v2f = v2(fl)
            core, quot, t, rms = abacus_data_4(lam)
            gc = G_of(core); vpc = vpi(*gc)
            fc = hook_f(core); v2fc = v2(fc)
            res = None
            if vpc != float('inf') and vp != float('inf'):
                res = vp - v2f - vpc + v2fc
            rows.append(dict(m=m, lam=lam, re=re, im=im, vp=vp, v2f=v2f,
                             core=core, quot=quot, t=t, rms=rms,
                             vpc=vpc, v2fc=v2fc, res=res,
                             qsz=sum(sum(q) for q in quot), csz=sum(core)))
    return rows

def unordered_quot(quot):
    return tuple(sorted(quot, key=lambda p: (sum(p), p)))

def group_test(rows, keyfn, label):
    g = defaultdict(set)
    detail = defaultdict(list)
    for r in rows:
        if r['res'] is None:
            continue
        g[keyfn(r)].add(r['res'])
        detail[keyfn(r)].append((r['lam'], r['res']))
    nc = [(k, sorted(v)) for k, v in g.items() if len(v) > 1]
    print(f"  [{label}]  #groups={len(g)}  non-constant={len(nc)}")
    return g, nc, detail

def main():
    MMAX = 11
    print(f"Building residual table for all lambda |- 2m, m<= {MMAX} ...", flush=True)
    rows = build_rows(MMAX)
    print(f"  total shapes: {len(rows)}")

    # ---- cross-checks ----
    oksize = all(r['csz'] + 4 * r['qsz'] == sum(r['lam']) for r in rows)
    inj = {}
    okinj = True
    for r in rows:
        key = (r['core'], r['quot'])
        if key in inj and inj[key] != r['lam']:
            okinj = False
        inj[key] = r['lam']
    van = [r['lam'] for r in rows if r['vp'] == float('inf')]
    # runner-multiset consistency with quotient size
    okrms = all(len(r['rms']) == r['qsz'] for r in rows)
    print(f"  cross-checks: core+4|quot|=n {oksize}; (core,quot) injective {okinj}; "
          f"runner-ms size matches {okrms}; full vanishers {van}")

    print("\n" + "=" * 78)
    print("LADDER of grouping keys (does runner data make residual a function?)")
    print("=" * 78)
    g0, nc0, _ = group_test(rows, lambda r: (r['vpc'], r['v2fc'], unordered_quot(r['quot'])),
                            "K0  (vpc, v2fc, UNORDERED quotient)")
    g1, nc1, _ = group_test(rows, lambda r: (r['vpc'], r['v2fc'], r['quot']),
                            "K1  (vpc, v2fc, ORDERED quotient)  [= K0 + which-runner]")
    g2, nc2, d2 = group_test(rows, lambda r: (r['vpc'], r['v2fc'], r['quot'], r['rms']),
                             "K2  (vpc, v2fc, ORDERED quotient, runner-multiset)")
    g3, nc3, d3 = group_test(rows, lambda r: (r['vpc'], r['v2fc'], r['quot'], r['t']),
                             "K3  (vpc, v2fc, ORDERED quotient, bead-counts t_j)")
    g4, nc4, _ = group_test(rows, lambda r: (r['core'], r['quot']),
                            "K4  (core, ORDERED quotient)  [complete invariant; sanity]")

    print("\n  --> runner effect (K0 -> K1): non-constant", len(nc0), "->", len(nc1),
          f"  resolved {len(nc0)-len(nc1)} of {len(nc0)}")
    print("  --> runner-multiset (K1 -> K2): non-constant", len(nc1), "->", len(nc2))
    print("  --> bead-counts t_j (K1 -> K3): non-constant", len(nc1), "->", len(nc3))

    print("\n" + "=" * 78)
    print("EXCEPTIONS that survive K1 (same vpc,v2fc,ORDERED quotient; diff residual)")
    print("=" * 78)
    g1d = defaultdict(list)
    for r in rows:
        if r['res'] is None:
            continue
        g1d[(r['vpc'], r['v2fc'], r['quot'])].append(r)
    surv = 0
    for k, rs in g1d.items():
        resset = set(r['res'] for r in rs)
        if len(resset) > 1:
            surv += 1
            print(f"  key vpc={k[0]} v2fc={k[1]} quot={k[2]}  residuals={sorted(resset)}")
            for r in rs:
                print(f"      lam={r['lam']} core={r['core']} t_j={r['t']} "
                      f"rms={r['rms']} res={r['res']}")
            if surv >= 8:
                print("   ... (truncated)")
                break

    print("\n" + "=" * 78)
    print("ZOOM: core=(5,2,1,1,1), single-box quotient on each runner")
    print("=" * 78)
    target_core = (5, 2, 1, 1, 1)
    found = [r for r in rows if r['core'] == target_core and r['qsz'] == 1]
    for r in sorted(found, key=lambda r: r['rms']):
        runner = r['rms'][0] if r['rms'] else None
        print(f"  lam={r['lam']}  box on runner {runner}  t_j={r['t']}  "
              f"vp={r['vp']} v2f={r['v2f']} residual={r['res']}")

    # ---- if K1 (ordered quotient) already a function, conjecture Delta(runner) ----
    print("\n" + "=" * 78)
    print("CONJECTURED correction: residual at fixed core vs runner of a single box")
    print("=" * 78)
    # For every core with several single-box placements, tabulate residual by runner
    bycore = defaultdict(dict)
    for r in rows:
        if r['res'] is None or r['qsz'] != 1:
            continue
        runner = r['rms'][0]
        bycore[r['core']].setdefault(runner, set()).add(r['res'])
    # show a few cores where residual varies with runner
    varying = [(c, d) for c, d in bycore.items()
               if len(set(next(iter(v)) for v in d.values())) > 1]
    print(f"  cores with single-box residual varying by runner: {len(varying)}")
    for c, d in sorted(varying, key=lambda x: sum(x[0]))[:12]:
        row = {runner: sorted(v) for runner, v in sorted(d.items())}
        print(f"    core={c}: residual by runner {row}")

    # ---- CSV ----
    csv_path = os.path.join(RESULTS, "residual-by-runner.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["m", "lam", "core", "ordered_quotient", "runner_multiset",
                    "bead_counts_t", "vp", "v2f", "vpc", "v2fc", "qsz", "residual"])
        for r in rows:
            w.writerow([r['m'], r['lam'], r['core'], r['quot'], r['rms'], r['t'],
                        r['vp'], r['v2f'], r['vpc'], r['v2fc'], r['qsz'], r['res']])
    print(f"\n  wrote {csv_path}  ({len(rows)} rows)")
    return rows, (nc0, nc1, nc2, nc3, nc4)

if __name__ == "__main__":
    main()
