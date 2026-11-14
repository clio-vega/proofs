# LEAN 2026-11-14 c2: Q559 `cor:branch` (the w_n → ∞ branching rule) is sorry-free

**Target.** `cor:branch` of `proofs/2026-10-23-Q559-vanishing-converse.tex` L427–L437:
under the hypotheses of thm:main with max m ≥ 0,
lim_{w_n→∞} w_n^{−(m_{j*}−r+1)} Ψ^C_{(1^r),θ}(M) = (−1)^{r+j*} Ψ^C_{(1^{r−1}),θ̂}(M∖{w_n}).
Registry node `clio-Q559-branching-rule`.
**Project.** `lean/tworow_d4_kernel` (github.com/clio-vega/tworow-d4-kernel), new module
`TworowD4Kernel/CorBranch.lean`. Commit **`f6e2c6b`** is pushed, and the upstream SHA was confirmed by
`gh api`, with an HTTP 422 negative control on a fake SHA. Full `lake build` rc=0 (read from `lake`
itself), 3489 jobs (ThmMain's run had 3488, so the new module is built). CorBranch has **0 warnings** by
`awk`. The same `awk` caught the one CorBranch `show` linter warning in the previous build, so that zero
is a real measurement.

## Encoding (decided before tactics): an equality, with no filters
Substitute w_n = u⁻¹. The paper's coefficientwise limit in K((w_n⁻¹)) is "the value at u = 0 of the
regular part". With LemOrder's `RegAt0` (f = P/Q wherever Q ≠ 0, with Q(0) ≠ 0), there are two statements:

```
theorem cor_branch_exists (w) (M) (hn : n ∈ M) (k) (θ : Fin (k+1) → ℤ) (js)
    (hadm : Admissible w (M.erase n)) (hinj : Injective (rhoShiftPaper θ))
    (harg : ∀ j, m j ≤ m js) (hpos : 0 ≤ m js) :
    ∃ F, RegAt0 F ∧ F 0 = (-1)^(k+js) * psiHead w (M.erase n) k (thetaHat θ js) ∧
      ∀ u, Admissible (update w n u⁻¹) M → u^(m js - k) * psiHead (update w n u⁻¹) M (k+1) θ = F u

theorem cor_branch [Infinite K] ... (G) (hG : RegAt0 G) (hGΨ : same agreement with G) :
    G 0 = (-1)^(k+js) * psiHead w (M.erase n) k (thetaHat θ js)
```

The second statement is the *limit*. It says **every** regular function that agrees with
u^{m*−k}Ψ has the paper's value at 0, so the value is determined by Ψ and is not an artefact of a
chosen F. Uniqueness is `RegAt0.zero_eq_of_eq_off`: over an infinite field, two regular-at-0 functions
that agree off a finite set agree at 0. It is a 3-line consequence of `RegAt0.exists_ne`. I did not
use `Filter.Tendsto` because ℚ(w_1,…,w_{n−1}) has no topology in Mathlib, and the equality form
exists.

The sign: `js` is 0-indexed, so (−1)^{k+js} = (−1)^{(k+1)+(js+1)} = (−1)^{r+j*}.

**Deviation (more general, not weaker):** |M| ≥ r is not assumed, since the identity never uses it.
n is any element of M (`glueOrd`, as in ThmMain), and w only needs to be admissible on M∖{n}.
For the last point I added `hBlock_zero'`, which is `hBlock_zero` with the restricted hypothesis.
The original proof only ever used the restriction.

## What builds sorry-free (sorries: **0**)
- `RegAt0.zero_eq_of_eq_off`: well-definedness of the limit.
- `badSet` / `admissible_update_inv`: off {0} ∪ {1/w_i : i ∈ M'}, the point w[n ↦ u⁻¹] is admissible.
- `order_identity_update`: u^{m_js−k} Ψ(w[n↦u⁻¹]) = orderF w … u at every admissible u. This is
  `psiHead_order_identity` transported back to w by `hBlock_congr`/`gBlock_congr`. It needs neither
  distinctness nor |M| ≥ r.
- `limit_eq_orderF_zero`: any agreeing regular G has G 0 = H(0). It needs only argmax and m* ≥ 0, which
  is why the lesions can use it.
- `cor_branch_exists` and `cor_branch`.

No new mathematics was needed. The pointwise ingredients from LemOrder/BlockA were already the
whole proof; what was missing was the uniqueness lemma that turns "F(0)" into "the limit".

## Paper-convention instance (no existential)
M = {0,1}, n = 1, w_0 = 2, θ = (0,0), so m = (0,1) and js = 1 (the paper's j* = 2). The exponent is 0,
the sign +1, θ̂ = (0).
- `inst_branch_value`: `orderF wI {0,1} 1 1 ![0,0] 1 (snocOrd 1 ordI') ordI' 0 = 3/4`, computed directly.
- `inst_cor_branch`: every agreeing regular G has G 0 = 3/4.
- An **independent** SymPy check takes the `limit` W→∞ of the head-injection sum eq:Psidef (it does not
  use the alternant or Lean). It gives **3/4** at θ=(0,0) and **−3/2** for W⁻¹Ψ at θ=(2,0), which matches
  the lesion instances below. It also reproduces Ψ({2,3}) = 5/9 = `inst_value`.

## Lesions (Lean derives False). Every non-lesioned hypothesis is discharged
`CorBranchWith Pinj Parg sgn` is `cor_branch_exists` over ℚ, ι = ℕ. The hypotheses are weakened to
`P ∨ _`, and the value is scaled by `sgn`.
- `corBranchWith_control`: `CorBranchWith False False 1` **holds**. It is `cor_branch_exists`.
- **L-dist** (drop distinctness): **refuted** at θ = (1,0), m = (1,1), js = 0 (an argmax). The true
  limit is 0, because the tied Laplace terms cancel. The claim is −3/2.
- **L-sign** (flip (−1)^{k+js}): **refuted** at θ = (2,0). The true limit is −3/2 and the claim is +3/2.
- **L-arg** (js need not be an argmax): **refuted** at θ = (0,0), js = 0. No claimed value is even
  needed: a regular F with u⁻¹Ψ = F(u) would make u·F(u) regular and agreeing with the true u⁰Ψ, so
  0 = H_{js=1}(0) = 3/4. *The exponent alone is impossible.* This is the mechanism by which the
  argmax hypothesis is load-bearing. It differs from BlockA's L-arg (wrong value of H(0)): here the
  wrong exponent makes the limit not exist at all.
- **Pair** (dist ∧ arg): refuted a fortiori.
- **Not run:** lesion of m* ≥ 0. With m* < 0 the G-block term u^{m*+1}G(u) has exponent ≤ 0, so the
  true limit picks up a Ψ^- contribution (at r=1, θ=(−1), the hand value is 1 − 3/4 = 1/4 ≠ 1, and the SymPy limit of W·Ψ at w_0=2 confirms **1/4**).
  `psiHead_order_identity` takes `hpos` as a hypothesis, so refuting this lesion would need that
  identity re-proved without it. Recorded, not done.

## `#print axioms` (same run as a `sorryAx` control)
```
cor_branch, cor_branch_exists, RegAt0.zero_eq_of_eq_off, limit_eq_orderF_zero,
inst_branch_value, inst_cor_branch, corBranchWith_control,
lesion_{dist,sign,arg,pair}_refutable        : [propext, Classical.choice, Quot.sound]   (11/11)
sorryControl1114c2                           : [sorryAx]
```

## Registry (`proofs/registry/rick-q488-converse.json`, backup `.bak-1114c2-lean`)
- **`clio-Q559-branching-rule` is PROMOTED `proved` → `lean-verified`**, with `lean` set to
  `CorBranch.cor_branch, cor_branch_exists, RegAt0.zero_eq_of_eq_off`.
- `clio-Q559-thm-main-formal`'s scope text is updated: cor:branch is no longer listed as unformalised.
- **Top node `clio-Q553-vanishing-converse-Q559` is NOT promoted**, because the ℚ(w) ↔ pointwise bridge
  is open.
- `registry_validate --proofs-dir .` reports **OK**. Planted control (an invalid trust value plus a
  missing file on this node): **3 problems**. Both plants fired, and the 3rd line is the boundary rule
  cascading up.

## Scoping note only (brief "Not today"): the ℚ(w) ↔ pointwise bridge
Local Mathlib has `MvPolynomial.funext` (`Mathlib/Algebra/MvPolynomial/Funext.lean:84`). Over an infinite
integral domain, p = q if eval x p = eval x q for **all** x. There is also `funext_set` (L58), which only
needs agreement on a box Π s_i with every s_i infinite. Admissibility (w_i ≠ 0, pairwise distinct) is
**not** a box condition, so neither applies directly. The natural route is to clear denominators.
Write Ψ = P/D with D = (product of the eq:Psidef denominators) ∈ ℚ[w]. Then D·P vanishes at every
point: at admissible points because Ψ does, and at inadmissible ones because D does. So `funext` gives
D·P = 0, and therefore P = 0 since D ≠ 0 in a domain. The real cost is *defining* Ψ as an element of
`FractionRing (MvPolynomial ι ℚ)` and proving that `eval` commutes with it at admissible points. That
is a new, sizeable encoding layer, and it was not started.
