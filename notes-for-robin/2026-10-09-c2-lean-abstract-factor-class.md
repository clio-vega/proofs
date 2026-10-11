# The factor class is a hypothesis now — and it wanted nonvanishing, not positivity

**LEAN 2026-10-09 c2.** Repo `clio-vega/tworow-d4-kernel`, commits `b259e81` + `d4cd3fe`, both
pushed and verified ancestors of `origin/main`.

File: [`TworowD4Kernel/AbstractFactorObstruction.lean`](https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/AbstractFactorObstruction.lean)
Snapshot: [`NOTES-2026-10-09-abstract-factor-class.md`](https://github.com/clio-vega/tworow-d4-kernel/blob/main/NOTES-2026-10-09-abstract-factor-class.md)

## What it says

`TworowD4Kernel.ProductForm.not_exists_abstract_product_form` — for `3 ≤ b`, `Odd b`, there is no
`c` and no finite family `f : Multiset ((ℝ → ℝ) × ℤ)` with every `f_i` nonvanishing on `(−1,0)`
such that `t^b − t^{b−1} + 1 = t^c ∏_i (f_i t)^{e_i}` for all real `t`.

The three existing Lean nodes each widened the **exponents** (`+1`, then all of `ℤ`). All three
hard-code the factor shape `1 − t^d`. This one makes the factor class a **hypothesis**.

**12 declarations, all sorry-free, `[propext, Classical.choice, Quot.sound]` on all 12.**

## The thing worth your time

My own gap note predicted that **positivity at the root** is the load-bearing hypothesis — on the
good evidence that the `±1` sign hypothesis had turned out unused. Formalising shows that is one
notch too strong. The chain is `zpow_ne_zero`, then `Multiset.prod_ne_zero`, then `t ≠ 0`, and
**no step of it mentions an order.** What the mechanism consumes is *nonvanishing*. Positivity is
only how `1 − t^{d_i}` happens to achieve it on `(−1,0)`.

The gap between the two is not empty, and I made that a theorem rather than a remark:
`factor_class_strictly_wider` — `f t = t` is nonvanishing at **every** point of `(−1,0)` and
positive at **none** of it. The payoff is immediate: `Φ₁(t) = t − 1` is negative there, so
`not_product_form_cyclotomic_low` (no `t^c (t−1)^{e₁}(t+1)^{e₂}` form) is unreachable from a
positivity hypothesis and from both `1 − t^d` statements.

All three existing theorems are re-derived as instances (primed names) **and** left standing in
their own files, proved directly. Those are different facts and I wanted both on the record.

One honest negative: Mathlib already had the packaging I was told to check for —
`Multiset.prod_ne_zero` / `Multiset.prod_pos` are exactly "product of nonzeros/positives at a
point", and they were already the engines of the two earlier files. So the core lemma is a two-line
*instance*, not a reproof. The abstraction cost almost nothing, which is itself the finding: the
earlier files were more specific than their proofs required.

## An instrument fault you may care about, since it is in my standing brief

My brief says to treat `declaration uses 'sorry'` as mis-scoped. Measured today in one planted arm,
side by side:

- the pattern **as written in the brief**, straight quotes → **0**, with a sorry live. A constant
  function and a false green. Lean 4.30.0 emits the message with **backticks**.
- the backtick pattern → **1** (one plant site).
- `#print axioms` → **8 contaminated declarations**.

Both faults are live at once: wrong delimiter, and then still mis-scoped. The grade rests on
`#print axioms` alone. Planted at the dependency *root*, the propagation was 8 of 9 with
`factor_class_strictly_wider` **spared** — correctly, it is the one declaration not routing through
the core lemma. `lake build` exit **0** with the sorry live.

## Two scope lines I had to correct mid-session

My brief (15:06) said the even-`b` case "needs a coefficient argument with no root location in it"
and assigned it to the same day's PROVE slot. **That slot finished at 20:31 and closed it** — by a
unit-circle argument, not coefficients. And `Φ₆` is now known to be the only cyclotomic that can
divide `D_b`. So the general-cyclotomic Lean instance is an **unformalised known result** and a
concrete next target, not an open question. Both corrections are in the file, the snapshot and the
registry node.

Nothing here formalises Theorem D itself: `Y^λ_ρ`, `P_λ`, Kostka–Foulkes and charge have no Lean
definitions in this project, and `t^b − t^{b−1} + 1` is simply written down. The parent registry
grade stays `proved`.
