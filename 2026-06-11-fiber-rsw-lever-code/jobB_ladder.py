"""
JOB B (arms PROVE Route A) — the e_2-rung valuation ladder + anchor + SYT model.

  delta_j := v2(M_{j+1}) - v2(M_j)   (per-rung change as we multiply by e_2).
  M_j = <s_lambda, h1^{2(m-j)} e2^j>,  M_0 = f^lambda.

CODE.md asks for:
  1. results/e2-rung-ladder.csv: delta_j per rung, and the 8-step sum on pure-M toggle pairs.
     (NB: the 8-step sum = -4 on J*-PAIRS is TAUTOLOGICAL -- proved this cycle; we still tabulate
      delta_j, which is genuine data, and flag the tautology honestly.)
  2. Decompose delta_j into binomial-shadow vs M-intrinsic; report which carries structure.
  3. Confirm the anchor v2(M_0)=v2(f^lambda) via the 2-quotient hook-length factorisation.
  4. Confirm Tool C (the positive SYT-count model M_j = sum_mu K^{(j)}_{lambda mu} f^mu via
     vertical-2-strip chains) -- explains WHY no global v2(M_j) formula (v2 of a positive sum).
"""
import sys, math, csv, os
from collections import Counter, defaultdict
import symfunc as sf
from job1_tie_census import partitions, M_vector, v2, analyze
from core_quotient import core_quotient

RESULTS = "/home/clio/projects/code/results"

# -------------------------------------------------- f^lambda (SYT count) via hook length
def f_lambda(lam):
    lam = [x for x in lam if x > 0]
    n = sum(lam)
    conj = [sum(1 for r in lam if r > c) for c in range(lam[0])] if lam else []
    prod = 1
    for i, row in enumerate(lam):
        for j in range(row):
            arm = row - 1 - j; leg = conj[j] - 1 - i
            prod *= arm + leg + 1
    return math.factorial(n) // prod

# -------------------------------------------------- (1)+(2) delta_j ladder
def run_ladder(MMAX):
    print(f"=== (1)/(2) e_2-rung ladder delta_j  (m<={MMAX}) ===")
    rows = []
    # 8-step sum on pure-M pairs (tautological cross-check)
    pure8_ok = pure8_tot = 0
    delta_dist = Counter()
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            vM = [v2(x) for x in Ms]      # None where M_j=0
            r = analyze(lam, m)
            J = set(r['Jstar'])
            for j in range(m):
                if vM[j] is None or vM[j+1] is None:
                    dj = None
                else:
                    dj = vM[j+1] - vM[j]
                    delta_dist[dj] += 1
                # binomial-shadow part of the FULL val step: val(j+1)-val(j) = 1 + 2 dbin + 2 dM
                dbin = (v2(math.comb(m, j+1)) - v2(math.comb(m, j))
                        if (Ms[j] != 0 and Ms[j+1] != 0) else None)
                rows.append(dict(m=m, lam='|'.join(map(str, lam)), j=j,
                                 Mj=Ms[j], Mj1=Ms[j+1],
                                 v2Mj=(vM[j] if vM[j] is not None else 'inf'),
                                 v2Mj1=(vM[j+1] if vM[j+1] is not None else 'inf'),
                                 delta_j=(dj if dj is not None else ''),
                                 dbin=(dbin if dbin is not None else ''),
                                 in_Jstar_j=int(j in J), in_Jstar_j1=int((j+1) in J)))
            # 8-step sum on pure-M toggle pairs {j, j+8} both in J*
            for j in J:
                if (j+8) in J and (j ^ 8) == j+8:
                    pure8_tot += 1
                    if vM[j] is not None and vM[j+8] is not None:
                        if vM[j+8] - vM[j] == -4:
                            pure8_ok += 1
    with open(f"{RESULTS}/e2-rung-ladder.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"  wrote results/e2-rung-ladder.csv ({len(rows)} rows)")
    print(f"  delta_j distribution: {dict(sorted((k,v) for k,v in delta_dist.items()))}")
    print(f"  pure-M {{j,j+8}} J*-pairs with v2(M_j+8)-v2(M_j) = -4: {pure8_ok}/{pure8_tot}")
    print(f"    (TAUTOLOGICAL on J*-pairs -- val constant on J* forces it; see steplaw doc)")

# -------------------------------------------------- (3) anchor v2(M_0)=v2(f^lambda) via 2-quotient
def run_anchor(MMAX):
    print(f"\n=== (3) anchor: v2(M_0)=v2(f^lambda) and the 2-quotient factorisation  (m<={MMAX}) ===")
    ok_eq = tot = 0
    quot_ok = quot_tot = 0
    fail = []
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            fl = f_lambda(lam)
            tot += 1
            if Ms[0] == fl:
                ok_eq += 1
            else:
                if len(fail) < 5: fail.append((tuple(lam), Ms[0], fl))
            # 2-quotient factorisation: v2(f^lambda) = v2(multinomial(2m;2|q0|,2|q1|)) + v2(f^q0)+v2(f^q1)
            core, quot = core_quotient(lam, 2)
            if sum(core) == 0:    # empty 2-core
                q0, q1 = quot[0], quot[1]
                a0, a1 = sum(q0), sum(q1)
                multinom = math.factorial(2*m) // (math.factorial(2*a0) * math.factorial(2*a1))
                pred = v2(multinom) + (v2(f_lambda(q0)) or 0) + (v2(f_lambda(q1)) or 0)
                act = v2(fl)
                quot_tot += 1
                if (act or 0) == pred:
                    quot_ok += 1
    print(f"  M_0 == f^lambda:                              {ok_eq}/{tot}")
    if fail: print(f"    M_0 != f mismatches: {fail}")
    print(f"  v2(f^lambda) via 2-quotient factorisation:    {quot_ok}/{quot_tot} (empty 2-core shapes)")

# -------------------------------------------------- (4) Tool C: positive SYT-count model
def e2perp_chains_Mj(lam, m, j):
    """M_j via the e_2^perp model: M_j = sum_mu (#vertical-2-strip chains lambda->mu of length j) f^mu.
    e_2^perp s_nu = sum over mu with nu/mu a vertical 2-strip of s_mu.  Apply j times, then f^mu = <s_mu,h1^N>."""
    # represent current symmetric function as Counter over partitions (Schur basis), start s_lambda
    from collections import Counter as C
    cur = C({tuple(x for x in lam if x>0): 1})
    for _ in range(j):
        nxt = C()
        for nu, coef in cur.items():
            # vertical 2-strip removal: remove 2 boxes, no two in same row (i.e. a vertical strip),
            # equivalently nu/mu is a vertical 2-strip <=> mu obtained by removing 2 boxes from
            # distinct rows that leaves a partition.  Use e_2^perp = adjoint of mult by e_2.
            for mu in remove_vertical_2strip(nu):
                nxt[mu] += coef
        cur = nxt
    # f^mu count weighted
    total = 0
    for mu, coef in cur.items():
        total += coef * (f_lambda(mu) if sum(mu) > 0 else 1)
    return total

def remove_vertical_2strip(nu):
    """All mu obtained from nu by removing a vertical 2-strip (2 boxes in distinct rows, result a partition)."""
    nu = list(nu)
    L = len(nu)
    out = set()
    # choose two distinct rows i<k to remove one box each, keeping partition shape
    for i in range(L):
        for k in range(i+1, L):
            a = list(nu)
            a[i] -= 1; a[k] -= 1
            if all(a[r] >= a[r+1] for r in range(len(a)-1)) and all(x >= 0 for x in a):
                mu = tuple(x for x in a if x > 0)
                out.add(mu)
    return out

def run_syt_model(MMAX=7):
    print(f"\n=== (4) Tool C: M_j = sum_mu (#vert-2-strip chains) f^mu  (m<={MMAX}) ===")
    ok = tot = 0; fail = []
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            for j in range(m+1):
                pred = e2perp_chains_Mj(lam, m, j)
                tot += 1
                if pred == Ms[j]:
                    ok += 1
                elif len(fail) < 5:
                    fail.append((tuple(lam), j, pred, Ms[j]))
    print(f"  Tool C model matches M_j: {ok}/{tot}")
    if fail: print(f"    mismatches: {fail}")
    print(f"  => every M_j is a NONNEGATIVE sum of f^mu (no cancellation); v2 of a positive sum")
    print(f"     has no min-of-valuations formula -- this is WHY no global v2(M_j) closed form exists.")

if __name__ == "__main__":
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    run_ladder(MMAX)
    run_anchor(MMAX)
    run_syt_model(7)
