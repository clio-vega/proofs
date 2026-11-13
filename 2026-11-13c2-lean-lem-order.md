# LEAN 2026-11-13 c2: Q559 `lem:order` (= `prop:caseA`, the order comparison) is sorry-free

**Target.** `lem:order`, i.e. `prop:caseA` of `proofs/2026-10-23-Q559-vanishing-converse.tex` L284–L302 (the label
`lem:order` was added alongside `prop:caseA`, as the brief asked).
**Project.** `lean/tworow_d4_kernel` (github.com/clio-vega/tworow-d4-kernel), module `TworowD4Kernel/LemOrder.lean`.
Commits: **`c15c287`** (the content) and **`11d8675`** (wrapping long lines). Both are pushed, and the upstream SHA was confirmed by
`gh api`, with an HTTP 422 negative control on a fake SHA. `lake build` rc=0, 3487 jobs, 0 warnings in LemOrder.
For the warning count I used `awk`, which had found the 4 line-length warnings before the fix. The `grep -c` wrapper printed 0 on the same
`/tmp` log both before and after.

## Step 0: encoding decision (R2), recorded before any tactics
R2 is the elementary route. R1 would have retyped the whole `Bfactor` layer over `Frac ℚ[w]`. R2 needs only `Polynomial.roots`. The
predicate is

    RegAt0 f  :=  ∃ P Q : K[X], Q(0) ≠ 0 ∧ ∀ u, Q(u) ≠ 0 → f(u)·Q(u) = P(u).

Note that `u = 0` is **not** exempted, so the value f(0) is the value of the rational function at 0.
- **Closure lemmas** (`regAt0_*`): constants, `u`, `+`, `−`, `*`, `^n`, `/const`, finite `∑`/`∏`, and `/R(u)` for any polynomial
  with R(0) ≠ 0.
- **H and G are regular at 0, now as theorems** (`regAt0_hBlock`, `regAt0_gBlock`). Both block files had only "read [this] off the
  definition". That scope limit is now discharged.
- **The nonvanishing principle** (`RegAt0.exists_ne`): f(0) ≠ 0 and K infinite ⇒ f is nonzero off any finite set. The proof is
  P(0) = f(0)Q(0) ≠ 0, so P ≠ 0, and P and Q have finitely many roots.

**What "Ψ ≠ 0" means here.** The field K is infinite. If w is admissible on M and Ψ(M', θ̂)(w) ≠ 0, then moving w_n alone gives an
admissible point w₂ with Ψ(M, θ)(w₂) ≠ 0. For K = ℚ this says Ψ is not identically zero as a rational function. I chose the pointwise
form of the IH because an induction can consume it.

## What builds sorry-free
- `orderF` defines F(u) = H(u) + (−1)^{k+1} u^{(m*+1).toNat} G(u).
- `psiHead_order_identity` proves u^{m*−k} Ψ = F(u) at w_n = u⁻¹, from `psiHead_split` + `psiPlus_eq` + `psiMinus_eq`.
- `orderF_zero` proves F(0) = H(0). **This is the only place m* ≥ 0 enters.** The G-term carries u^{m*+1} with m*+1 ≥ 1.
- `lem_order` (main theorem): the hypotheses are `hadm`, `n ∈ M`, `OrdOK` on M and on M', `SnocCompat`, m injective, j* an argmax,
  m* ≥ 0 and the pointwise IH. The conclusion is ∃ v, Admissible (update w n v) M ∧ Ψ(update w n v) ≠ 0. The excluded finite set is
  {0} ∪ {1/w_i : i ∈ M'}. A finite congruence (`hBlock_congr`, `gBlock_congr`, `Bfactor_congr`) moves H and G from w₂ back to w.

Sorries: **0**.

## Lesions (Lean derives `False`). Every non-lesioned hypothesis is discharged on the instance
The instance is M = {0,1}, n = 1, k = 1, ord' ≡ 0 and ord = `snocOrd 1 ord'`. `OrdOK` is proved for this ord (`ordI_ok`), and
`SnocCompat` comes from `snocOrd_compat`. The statements are `LemOrderWith Pinj Pih`, which has the shape of `lem_order` with a
hypothesis weakened to `P ∨ h`.
- `lemOrderWith_control`: `LemOrderWith False False` **holds**. It is `lem_order` itself.
- **L-IH** (`LemOrderWith False True`): **refuted.** At w = (1,3), θ = (2,0) we have m = (2,1), which is distinct, with m* = 2. The
  diagonal factor 1 − w₀⁻² = 0 sits in the only T = {0,1}, so Ψ(M) = 0 for every w_n.
- **L-dist** (`LemOrderWith True False`): **refuted.** At w = (2,3), θ = (1,0) we have m = (1,1). **The IH holds there**
  (Ψ(M',θ̂) = 3/2), yet Ψ(M) ≡ 0 by the forward ρ-shift law.
- **Pair** (`LemOrderWith True True`): refuted a fortiori. Every pair contains a refutable single, so the joint-fatality question
  does not arise here.
- **L-gap** (drop m* ≥ 0): **refuted at the step it feeds** (`orderF_zero`). At M = {0,1}, k = 0, θ = (−1), m* = −1 we get
  F(0) − H(0) = −G(0) = −3/4. *As a statement*, `lem:order` without m* ≥ 0 is **not** refutable, because `thm:main` covers m* ≤ −1
  through the reflection. So the gap is load-bearing for this proof but not for the truth of the conclusion. Likewise, L-IH is
  refutable only because the lemma is pointwise. Globally the IH is supplied by induction.

## `#print axioms` (same run as a `sorryAx` control)
```
lem_order, psiHead_order_identity, RegAt0.exists_ne, regAt0_orderF, orderF_zero, lemOrderWith_control,
lesion_{IH,dist,pair,gap}_refutable                      : [propext, Classical.choice, Quot.sound]   (10/10)
sorryControl                                             : [sorryAx]
```

## Registry
I added a new node, `clio-Q559-order-comparison` (`lean-verified`), under `clio-Q559-two-block-leading-order`. **The parent is not
promoted.** `registry_validate --proofs-dir .` reports `OK`. Its positive control is a copy with an invalid trust value and a missing
file planted, and the validator reports 3 problems on it.

## Still unformalised (for the parent and for thm:main)
1. **Case B** (max m ≤ −1). Reflect with `lem_reflect` at the point w⁻¹ and use θ∨ (m∨ = −m). This needs the IH for θ∨'s θ̂, which comes
   from the global induction and not from the same point. So Case B cannot be stated pointwise like Case A: it needs the *global*
   form `∃ w, Admissible w M ∧ Ψ(w) ≠ 0`.
2. **The induction on r.** It needs:
   - the base case Ψ(M, 0, ()) = 1;
   - `exists_admissible` over an infinite K, by Finset induction;
   - a `psiHead` congruence, so that Ψ(M') does not see w_n;
   - an `OrdOK M (k+1)` ordering that satisfies `SnocCompat`. `snocOrd` alone is OrdOK only on T ∋ n, so the ordering has to be
     glued with any ordering on the T that avoid n;
   - injectivity of θ̂'s shifts (m ∘ succAbove).

   These are about 150 lines and contain no new mathematics. Next LEAN target: steps 1 and 2 together, which would complete thm:main.
