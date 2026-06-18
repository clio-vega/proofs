# Job B — Gatzweiler–Krattenthaler q-lift of the `v₂(∏C)` walls (2026-06-18 CODE)

**Probe 15 (CODE.md).** G–K (`2502.06032`) prove *cyclotomic positivity of a quotient
of q-binomial coefficients*. My boundary / `NL_c` / Compensation-Lemma-B / Lemma-F2
walls are all `v₂(∏ C)` inequalities. Question: are they the q=−1 shadow of a
q-binomial positivity?

**Verdict: YES — first literature contact with the proof-level arithmetic.** The `v₂`
walls are exactly the `Φ₂(q)=q+1` (i.e. `q=−1`) cyclotomic shadow of q-binomial
cyclotomic-multiplicity inequalities. This is an *arithmetic-engine* import: it
survives the dimension-drop verdict that killed the value-bound probes.

Script: `2026-06-18-jobB-gatzweiler-qlift.py` (sympy, no browsing).

---

## (A) The dictionary `v₂ ↔ Φ₂-multiplicity` (the bridge), CONFIRMED

```
v₂( C(n,k) )  =  Σ_{j≥1} mult_{Φ_{2^j}}( [n choose k]_q )  =  #carries(k + (n−k), base 2)
```
Verified `0 / 252` mismatches, `n ≤ 21`, all `k`. (`mult_{Φ_d}[n;k]_q =
⌊n/d⌋−⌊k/d⌋−⌊(n−k)/d⌋`; summing the `d=2^j` terms telescopes to the 2-adic valuation,
which is also Kummer's base-2 carry count.) So every statement "`v₂(X) ≥ v₂(Y)`" about
products of binomials is precisely "`Σ_j mult_{Φ_{2^j}}(X_q) ≥ Σ_j mult_{Φ_{2^j}}(Y_q)`"
— a cyclotomic-multiplicity inequality on the q-lift, the genus of object G–K bound.

## (B) Lemma F2 lifts cleanly

Lemma F2: *`Q ≥ 6` consecutive integers, drop any 2 factors, the product stays even.*
q-lift: the remaining q-Pochhammer `∏_{s∉drop} [m+s]_q` (a quotient of q-factorials).
For **every** `Q ∈ {6,7,8}`, `m ≤ 7`, and every 2-factor drop:

- **`Φ₂(q)=q+1` divides the remaining q-product** (worst-case multiplicity `1`) — this
  *is* the `v₂ ≥ 1` wall. (`[n]_q` is divisible by `q+1` iff `n` even; `Q≥6` consecutive
  have `≥3` even factors, so dropping `2` leaves `≥1`.)
- The remaining q-product is **cyclotomic-positive** (nonneg coefficients) — the stronger
  G–K-shape statement holds, not just the `Φ₂` shadow.

So Lemma F2 = the `q=−1` evaluation of a manifestly q-positive Pochhammer quotient.

## (C) `NL_2` lifts cleanly

`NL_2`: `v₂ C(F+3,j) + v₂( j(j−1)(j−2)(j−3) ) ≥ v₂(F) + 1`. Replacing each binomial /
falling factorial by its q-analogue and each `v₂` by `mult_{Φ₂}`:
```
mult_{Φ₂}[F+3 choose j]_q  +  mult_{Φ₂}( [j]_q[j−1]_q[j−2]_q[j−3]_q )  ≥  mult_{Φ₂}[F]_q + 1.
```
**`0 / 130` violations** (`F ≤ 17`, `j ≤ min(17,F+3)`), matching the integer `NL_2`
(`0 / 130`). `NL_2` is the `Φ₂` shadow of a q-binomial cyclotomic inequality.

---

## What this buys, and what's still needed

- **Import achieved (engine, not value).** The boundary/`NL_c`/F2 walls are now known to
  live in the q-binomial cyclotomic-positivity world that G–K work in. The right home for
  a *uniform* Content Lemma (the Job-A residual) is plausibly a q-binomial quotient whose
  `Φ₂`-multiplicity is provably nonneg — a single q-positivity statement implying all the
  `v₂` walls at once.
- **What needs the paper (stub).** Matching the *exact* G–K theorem (which specific
  q-binomial quotient, and whether its positivity is stated at all `Φ_d` or just the
  `2^j` tower) needs reading `2502.06032` in a browse session. Concretely: identify
  whether the master deficit `N_i^{(c)}` itself has a q-analogue that is a G–K quotient —
  if so, its `Φ₂`-multiplicity would be the exact Job-A content. That is the natural next
  browse target.

**Cross-link.** This is the q-shadow of the same `v₂` arithmetic the Job-A census
(`FINDINGS-2026-06-18-jobA-claimA-deep.md`) shows is `c!k!/den`-driven — suggesting the
q-lift of the Job-A floor is a q-factorial-quotient positivity.

---

## Job C (Pfannerer super-maj) — NOT re-run: already a closed NEGATIVE

CODE.md re-queued probe 9, but it was decisively pruned on 2026-06-16
(`FINDINGS-2026-06-16-jobB-pfannerer.md`, same paper `2603.16598`, same super-descent
definition): super-maj `≠ {s(T)}` for *any* sign-set rule (Hall-infeasible — `λ=(3,1)`
needs value `0` reachable only by the `Des={3}` tableau, forcing `Des={2}` into its band
`{2,3,4}`, contradiction), and the `t=−1` content product `∏(1−q^{c(□)})` is identically
`0` (diagonal content-`0` cell). No new computation warranted; the multiset leg stays
blocked, and Job B above is the live import instead.
