"""
JOB A (2026-06-11) -- Albion's z-asymmetric 4-core/4-quotient vs the d=4 residual.

Albion 2501.18520 ("Character factorisations, z-asymmetric partitions and
plethysm"): the Littlewood decomposition lambda <-> (t-core, t-quotient) via the
t-Maya diagram. A t-CORE mu is encoded by its CHARGE VECTOR
    kappa_t(mu) = (c_0, ..., c_{t-1}),  c_r = (#beads on runner r) - N/t,  sum c_r = 0,
which is the runner-imbalance.  z-asymmetric partitions lambda=(u+z|u) are exactly
those whose charge + quotient obey the conjugacy relations c_r + c_{z-r-1}=0 etc.
(Thm 2.3).  The "z-asymmetric quotient datum" the CODE plan asks about is therefore
    (charge vector kappa_4(4-core),  ordered 4-quotient (lam^0,...,lam^3)).

GOAL.  Does that datum supply a CLOSED FORM for
    residual(lambda) = v_pi(G_lambda) - v2(f^lambda) - v_pi(G_core) + v2(f_core)?
Decisive sub-questions:
  (Q0) Is kappa_4 a complete invariant of the 4-core? (cores <-> charge bijection)
       -> if yes, (kappa_4, quotient) <-> lambda is the COMPLETE INVARIANT, so
          "residual = f(kappa_4, quotient)" is TRUE-BUT-VACUOUS (the K4 tautology).
  (Q1) Does it SEPARATE the n=10 counterexample (1^6) vs (4,4,2)?  (charge differs?)
  (Q2) GENUINE content: is residual a function of a STRICTLY SMALLER datum --
        (charge, quotient SIZES only)?  (charge, multiset of quotient parts)?
        -- i.e. does the charge vector let us forget which runner / the box layout?
  (Q3) Is residual ADDITIVE / closed-form over the charge+quotient?
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
            prod *= (row - (j + 1)) + (conj[j] - (i + 1)) + 1
    return factorial(n) // prod

def beta_set(lam, r):
    lam = list(lam) + [0] * (r - len(lam))
    return [lam[i] + (r - 1 - i) for i in range(r)]

def beta_to_part(betas):
    betas = sorted(betas, reverse=True)
    r = len(betas)
    part = [betas[i] - (r - 1 - i) for i in range(r)]
    return tuple(x for x in part if x > 0)

def abacus4(lam, N=None):
    """4-core, ordered 4-quotient, canonical CHARGE vector kappa (sum 0)."""
    lam = [x for x in lam if x > 0]
    if N is None:
        N = len(lam)
        while N % 4 != 0:
            N += 1
        N = max(N, 4) + 8            # generous, multiple-of-4 headroom
        N -= N % 4
    betas = beta_set(lam, N)
    runners = {j: [] for j in range(4)}
    for b in betas:
        runners[b % 4].append(b // 4)
    quotient = tuple(beta_to_part(runners[j]) for j in range(4))
    t = tuple(len(runners[j]) for j in range(4))          # bead counts
    mean = N // 4
    charge = tuple(t[j] - mean for j in range(4))          # canonical, sum 0
    core_betas = [4 * k + j for j in range(4) for k in range(t[j])]
    core = beta_to_part(core_betas)
    return core, quotient, charge

@lru_cache(maxsize=None)
def Gtab(m):
    return psi_power_schur(m)

def G_of(part):
    n = sum(part)
    if n == 0:
        return (1, 0)
    return Gtab(n // 2).get(tuple(part), (0, 0))

def build(MMAX):
    rows = []
    for m in range(1, MMAX + 1):
        for lam, (re, im) in Gtab(m).items():
            vp = vpi(re, im)
            v2f = v2(hook_f(lam))
            core, quot, charge = abacus4(lam)
            gc = G_of(core); vpc = vpi(*gc); v2fc = v2(hook_f(core))
            res = None
            if vpc != float('inf') and vp != float('inf'):
                res = vp - v2f - vpc + v2fc
            rows.append(dict(m=m, lam=lam, vp=vp, v2f=v2f, core=core, quot=quot,
                             charge=charge, vpc=vpc, v2fc=v2fc, res=res,
                             qsz=sum(map(sum, quot))))
    return rows

def nonconst(rows, keyfn):
    g = defaultdict(set)
    for r in rows:
        if r['res'] is not None:
            g[keyfn(r)].add(r['res'])
    nc = {k: sorted(v) for k, v in g.items() if len(v) > 1}
    return len(g), nc

def qsizes(quot):  return tuple(sum(q) for q in quot)
def qparts(quot):  return tuple(tuple(sorted(q, reverse=True)) for q in quot)
def qmultiset(quot): return tuple(sorted(qsizes(quot)))

def main():
    MMAX = 9                       # n <= 18 as the plan asks
    print(f"Building residual table, all lambda |- 2m, m<={MMAX} (n<={2*MMAX}) ...",
          flush=True)
    rows = build(MMAX)
    used = [r for r in rows if r['res'] is not None]
    print(f"  shapes={len(rows)}  with finite residual={len(used)}\n")

    # ---------- Q0: charge vector <-> 4-core bijection? ----------
    print("=" * 74)
    print("Q0  Is the charge vector a COMPLETE invariant of the 4-core?")
    print("=" * 74)
    core2charge = {}
    charge2core = {}
    bij = True
    for r in rows:
        c, ch = r['core'], r['charge']
        if c in core2charge and core2charge[c] != ch: bij = False
        if ch in charge2core and charge2core[ch] != c: bij = False
        core2charge[c] = ch; charge2core[ch] = c
    print(f"  distinct 4-cores seen: {len(core2charge)}")
    print(f"  core -> charge and charge -> core both single-valued: {bij}")
    print(f"  => (charge, quotient) is the COMPLETE invariant of lambda: {bij}")

    # ---------- Q1: separation of the n=10 counterexample ----------
    print("\n" + "=" * 74)
    print("Q1  Does the charge vector separate (1^6) vs (4,4,2)?")
    print("=" * 74)
    for r in rows:
        if r['lam'] in ((1,1,1,1,1,1), (4,4,2)):
            print(f"  lam={r['lam']}  core={r['core']}  charge={r['charge']}  "
                  f"quot_sizes={qsizes(r['quot'])}  residual={r['res']}")

    # ---------- Q2: ladder of data, coarse -> fine ----------
    print("\n" + "=" * 74)
    print("Q2  Which datum makes residual a FUNCTION?  (#non-constant groups)")
    print("=" * 74)
    tests = [
        ("ordered quotient ALONE", lambda r: r['quot']),
        ("charge ALONE", lambda r: r['charge']),
        ("(charge, quotient SIZES)", lambda r: (r['charge'], qsizes(r['quot']))),
        ("(charge, quotient size MULTISET)", lambda r: (r['charge'], qmultiset(r['quot']))),
        ("(charge, quotient PARTS per runner)", lambda r: (r['charge'], qparts(r['quot']))),
        ("(charge, ordered quotient)  [= complete inv]", lambda r: (r['charge'], r['quot'])),
        ("(4-core, ordered quotient)  [complete inv sanity]", lambda r: (r['core'], r['quot'])),
    ]
    for label, fn in tests:
        ng, nc = nonconst(used, fn)
        print(f"  {label:48s} groups={ng:5d}  non-constant={len(nc)}")

    # ---------- Q3: additivity over quotient boxes at fixed charge ----------
    print("\n" + "=" * 74)
    print("Q3  Is residual additive over the quotient at fixed charge?")
    print("    test: residual(charge,q) =? base(charge) + sum_r w(charge,r)*|q^r|")
    print("=" * 74)
    # base(charge): residual of empty quotient (lambda = core)
    base = {}
    for r in used:
        if r['qsz'] == 0:
            base[r['charge']] = r['res']
    # single-box increments w(charge, runner)
    w = {}
    consistent_single = True
    for r in used:
        if r['qsz'] == 1:
            runner = next(j for j in range(4) if sum(r['quot'][j]) == 1)
            key = (r['charge'], runner)
            inc = r['res'] - base.get(r['charge'], 0)
            if key in w and w[key] != inc:
                consistent_single = False
            w[key] = inc
    print(f"  single-box increment w(charge,runner) well-defined: {consistent_single}")
    # now test the additive prediction on ALL shapes
    add_ok = add_bad = 0
    bad_examples = []
    for r in used:
        ch = r['charge']
        if ch not in base:
            continue
        pred = base[ch]
        ok = True
        for j in range(4):
            sz = sum(r['quot'][j])
            if sz == 0:
                continue
            if (ch, j) not in w:
                ok = False; break
            pred += w[(ch, j)] * sz
        if not ok:
            continue
        if pred == r['res']:
            add_ok += 1
        else:
            add_bad += 1
            if len(bad_examples) < 12:
                bad_examples.append((r['lam'], ch, qsizes(r['quot']), pred, r['res']))
    print(f"  additive prediction base+sum w*|q^r|:  ok={add_ok}  bad={add_bad}")
    if bad_examples:
        print("  first failures (lam, charge, qsizes, pred, actual):")
        for e in bad_examples:
            print("    ", e)

    # ---------- write CSV ----------
    path = os.path.join(RESULTS, "jobA-albion-charge.csv")
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["m","lam","core","charge","quotient","quot_sizes","residual"])
        for r in rows:
            wr.writerow([r['m'], r['lam'], r['core'], r['charge'], r['quot'],
                         qsizes(r['quot']), r['res']])
    print(f"\n  wrote {path}")
    return rows

if __name__ == "__main__":
    main()
