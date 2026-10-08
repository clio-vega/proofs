# LEAN 2026-10-08 c2 — `not_signed_product_form`: the denominator half of Theorem D

**Target:** extend the Theorem D product-form obstruction from numerator-only exponents to
exponents of sign `±1`.
**Project:** `lean/tworow_d4_kernel` (Lean 4.30.0, Mathlib v4.30.0).
**New file:** `TworowD4Kernel/SignedProductFormObstruction.lean`, imported from the root.
**Status: closed, sorry-free.** 9 new declarations, root `lake build` exit 0 (3205 jobs).

---

## 1. Target declarations

| declaration | statement |
|---|---|
| `TworowD4Kernel.ProductForm.not_signed_product_form` | for odd `b ≥ 3`, no `(c : ℕ, d : Multiset (ℕ × ℤ))` with all `d_i ≥ 1` has `t^b − t^{b−1} + 1 = t^c ∏_i (1 − t^{d_i})^{e_i}` for all real `t`, with `e_i` ranging over **all of `ℤ`** |
| `not_product_form_pm_one` | the literal paper display, `e_i ∈ {+1, −1}` |
| `not_product_form_ratio` | the quotient form `t^c ∏(1−t^{n_i}) / ∏(1−t^{m_j})` |

Supporting: `one_sub_zpow_pos`, `prod_signed_pos`, `prod_signed_ne_zero`.
Controls: `cleared_lhs_eq_zero_at_one`, `signed_prod_eq_zero_of_exponent_zero`,
`not_signed_prod_ne_zero_without_hd`.

Citation carried into the file: Theorem D, `proofs/2026-10-07-two-part-green-polynomials.tex`;
the root input is `TworowD4Kernel.exists_root_Ioo` from `NonCyclotomicRoot.lean`.

## 2. The finding: the prescribed route was unnecessary **and its premise was false**

The session brief prescribed: *clear the denominators first*, reducing to
`D_b · ∏_{e_i=−1}(1−t^{d_i}) = t^c ∏_{e_i=+1}(1−t^{d_i})`, then separate the two sides by their
value at `t = 1` — "the LHS has nonzero value at `t=1` (`1 − 1 + 1 = 1` times a product of zeros —
check this carefully, it may be the whole proof or may collapse)".

It collapses, and not in the direction the parenthesis anticipated. `D_b(1) = 1`, but the *cleared*
left side is `D_b(1) · ∏_{k ∈ den}(1 − 1^k) = 1 · 0 = 0` whenever `den` is nonempty — which is
precisely the case the denominators were introduced to handle. The "product of zeros" makes the LHS
**zero**, not nonzero. So the proposed `t = 1` separation has no content exactly where it is needed.

That dead route is now a theorem rather than a memory, so it cannot be re-proposed:

```lean
theorem cleared_lhs_eq_zero_at_one (b : ℕ) (den : Multiset ℕ) (hden : den ≠ 0) :
    ((1:ℝ)^b - (1:ℝ)^(b-1) + 1) * (den.map (fun k => 1 - (1:ℝ)^k)).prod = 0
```

stated for **every** `b` and every nonempty `den` — no oddness, no `b ≥ 3`.

**What actually closes it** is the same root in `(-1,0)`, used once more. Nonvanishing is preserved
under inversion:

> At `t₀ ∈ (-1,0)` every factor `1 − t₀^{d_i}` is strictly **positive** (`one_sub_pow_pos`, reused
> from the `+1` file), and a positive real raised to any integer power is positive (`zpow_pos`). So
> each `(1 − t₀^{d_i})^{e_i} > 0`, the product is positive (`Multiset.prod_pos`), and `t₀^c ≠ 0`.
> The right side is nonzero at `t₀`; the left side is `0`.

No denominator clearing, no `t = 1`, no cyclotomic theory. `one_sub_zpow_pos` is a one-line
composition. This is the **second** time on this file that a gap note's named route was not the
route the step needed — yesterday it was `isRoot_cyclotomic_iff`. A gap note records the route its
author could see.

## 3. Generality is free; the `±1` hypothesis is not load-bearing

Positivity of the base is the only thing the argument consumes, and it is indifferent to the
exponent. So `ℤ` costs exactly what `±1` costs. Consequently `not_product_form_pm_one`'s hypothesis
`hpm : ∀ p ∈ d, p.2 = 1 ∨ p.2 = -1` is **literally unused** — bound as `_hpm`, and its docstring
says so. Deleting `hpm` gives `not_signed_product_form` verbatim, still true. The declaration exists
only so a reader matching the formalisation against the paper's display finds the display.

## 4. `#print axioms` — the grade

`grep "declaration uses 'sorry'"` is a **dead instrument** in this repo (it has read 0 with a sorry
planted), so the grade rests on `#print axioms`, run as a two-arm canary.

**Clean arm** — all 9 declarations:

```
'TworowD4Kernel.ProductForm.one_sub_zpow_pos'                depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.prod_signed_pos'                 depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.prod_signed_ne_zero'             depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.not_signed_product_form'         depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.not_product_form_pm_one'         depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.not_product_form_ratio'          depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.cleared_lhs_eq_zero_at_one'      depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.signed_prod_eq_zero_of_exponent_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.not_signed_prod_ne_zero_without_hd'   depends on axioms: [propext, Classical.choice, Quot.sound]
```

Standard three only. No `sorryAx`, no `native_decide`.

**Planted arm** — `sorry` substituted for the body of `one_sub_zpow_pos` (pattern asserted present,
predicted 1 / found 1, *before* mutating): **5** declarations carry `sorryAx`, along exactly the
chain `one_sub_zpow_pos → prod_signed_pos → prod_signed_ne_zero → not_signed_product_form →
not_product_form_pm_one`. The reading is **0 → 5**, so the instrument is alive; equal readings would
have meant it was dead. File restored, `diff -q` clean, root rebuilt exit 0, axioms re-measured at 0
`sorryAx`.

**Unplanned corroboration.** The planted sorry did *not* reach `not_product_form_ratio`. That
independently confirms that declaration's prose claim to be proved *directly* rather than
transported along `not_signed_product_form` — a propagation pattern measures the dependency graph,
not only cleanliness. I had asserted this in a docstring; the canary checked it for free.

**Comment-stripped grep.** Block comments removed with `awk` (12 docstrings in this repo say
"sorry" in prose): **0** occurrences of `sorry` outside comments.

## 5. Ablation — as a counterexample, not a build failure

`not_signed_prod_ne_zero_without_hd` proves that `prod_signed_ne_zero` with the `d_i ≥ 1` hypothesis
**deleted is false**: witness `d = {(0,−1)}` at `t = −1/2`, since `(1 − t^0)^{−1} = 0^{−1} = 0` in
Lean. A failed `lake build` would be a fact about my proof; this is a fact about the statement.

Note the *new* content over the `+1` file's `prod_eq_zero_of_exponent_zero`: a zero exponent now
trivialises the identity from the **denominator** side as well, because Lean's `0^{−1} = 0` and
`x / 0 = 0`. That convention is Lean's, not the paper's, and it is the reason `1 ≤ d_i` is doing
real work in every signed statement here. The file says so explicitly.

## 6. Overlap guard — statement shape, not name

Before writing anything I grepped the repo for the *statement shape*, not the name: `zpow` / integer
exponents (3 hits, all `QuantumInteger.lean`, unrelated), `Ioo (-1` (only `ProductFormObstruction`
and `NonCyclotomicRoot`), quotient-of-products shapes (no hits), and `product form` (the two other
hits, `TypedBlocks`/`LabelledBlocks`, are multinomial block counts, a different sense of the phrase).
No prior formalisation of the signed case existed.

## 7. Scope — what is still NOT formalised

Unchanged from the `+1` file, and nothing here should be read otherwise. **Theorem D itself is not
formalised.** `Y^λ_ρ`, Hall–Littlewood `P_λ`, Kostka–Foulkes polynomials and the charge statistic
have no Lean definitions in this project, and the identification
`D_{a,b} = Y^{(a,b)}_{(a,b)} = t^b − t^{b−1} + 1` is paper-side (Theorem C at `m = b`). Here
`t^b − t^{b−1} + 1` is simply *written down*. What this session adds is exponent generality in the
obstruction **mechanism**.

`unproved ≠ unformalised`. The parent registry node `thm-D-product-form-obstruction` stays
`proved`, not `lean-verified`.

## 8. Bookkeeping

- Registry: new node `lean-no-signed-product-form` under `thm-D-product-form-obstruction`,
  `trust: lean-verified`, `lean: TworowD4Kernel.ProductForm.not_signed_product_form`.
  `registry_validate.py --proofs-dir .` → 2 problems, **both pre-existing and on a different node**
  (`rick-two-point-formula-thm25`: `peer-claimed` trust, missing `rick/` sub-registry); 0 on mine.
  Without `--proofs-dir .` the validator double-prefixes and reports ~40 phantom missing files —
  that is the known flag hazard, not a regression.
- `ProductFormObstruction.lean`'s scope paragraph claimed the `±1` case was "*not* formalised
  below". That claim is now false, so it was amended to point at the new file (1/1 substitution,
  asserted before writing). Leaving it would have left a false statement in the repo.
