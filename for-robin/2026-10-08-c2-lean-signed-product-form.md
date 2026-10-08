# LEAN 2026-10-08 c2 — Theorem D's obstruction now covers denominators, and the route I was given was false

**One line:** `TworowD4Kernel.ProductForm.not_signed_product_form` is sorry-free — for odd `b ≥ 3`
there is no identity `t^b − t^{b−1} + 1 = t^c ∏_i (1 − t^{d_i})^{e_i}` with `d_i ≥ 1` and `e_i`
ranging over **all of ℤ**, not merely `±1`. That closes the scope gap yesterday's numerator-only
result declared in its own docstring.

- Lean: https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/SignedProductFormObstruction.lean
- Report: https://github.com/clio-vega/proofs/blob/main/2026-10-08-lean-not-signed-product-form.md

9 declarations, axioms `[propext, Classical.choice, Quot.sound]`, no `sorryAx`, no `native_decide`,
root `lake build` exit 0.

## The thing worth your time

The session plan told me to clear the denominators and then separate the two sides by their value at
`t = 1` — "the LHS has nonzero value at `t=1` (`1 − 1 + 1 = 1` times a product of zeros — check this
carefully, it may be the whole proof or may collapse)."

**The premise is false, and the sentence stating it contains its own refutation.** "1 times a product
of zeros" *is zero*. `D_b(1) = 1`, but the cleared left side is `D_b(1) · ∏_{k∈den}(1 − 1^k) = 0` as
soon as the denominator set is nonempty — which is exactly the case denominators were introduced to
handle. So the proposed separation has no content precisely where it is needed.

I have a standing lesson from yesterday that *a stated route is a claim about necessity*. This is one
rung further: the route was not merely unnecessary, its **premise was false**, and it was refutable
by arithmetic visible in the half-sentence that proposed it. I did not need to compute anything to
find that — I needed to read the clause instead of acting on its conclusion.

It is now a theorem, so it cannot be re-proposed from memory:

```lean
theorem cleared_lhs_eq_zero_at_one (b : ℕ) (den : Multiset ℕ) (hden : den ≠ 0) :
    ((1:ℝ)^b - (1:ℝ)^(b-1) + 1) * (den.map (fun k => 1 - (1:ℝ)^k)).prod = 0
```

Every `b`, every nonempty `den` — no oddness, no `b ≥ 3`.

## What actually closes it (four lines, as predicted — just not those four)

The same root in `(-1,0)` does it again, because **nonvanishing is preserved under inversion**. At
`t₀ ∈ (-1,0)` each factor `1 − t₀^{d_i}` is strictly *positive* (`one_sub_pow_pos`, reused verbatim
from yesterday), and a positive real to any integer power is positive (`zpow_pos`). Product positive,
monomial nonzero, so the right side is nonzero at `t₀` while the left side is `0`.

No denominator clearing, no `t = 1`, no cyclotomic theory — the same shape as yesterday's finding
that `isRoot_cyclotomic_iff` was a detour. Positivity of the *base* is all the argument consumes and
it is indifferent to the exponent, which is why `ℤ` came free and why `not_product_form_pm_one`'s
`±1` hypothesis is **literally unused** (bound `_hpm`, flagged in its docstring). I kept that
declaration only so a reader matching Lean against your display finds the display.

## Two small instrument notes

**The canary measured the dependency graph, not just cleanliness.** Planting a `sorry` in the most
upstream lemma took `#print axioms` from 0 to 5 `sorryAx` — and it did *not* reach
`not_product_form_ratio`. I had claimed in that declaration's docstring that it was proved directly
rather than transported along `not_signed_product_form`. The canary confirmed the claim for free. A
propagation *pattern* carries more information than a propagation *count*.

**`git status` on the parent repo read clean while my new file sat untracked inside it.**
`lean/.gitignore` line 6 is `tworow_d4_kernel/`, because that subdirectory is its own repo. The
parent's "clean" was true and about the wrong repo; pushing `git -C .../lean` would have been a
green no-op. Caught by asking why a repo I had just written to was clean, rather than reading the
word.

## Scope, unchanged

**Theorem D itself is still not formalised.** `Y^λ_ρ`, Hall–Littlewood `P_λ`, Kostka–Foulkes and
charge have no Lean definitions in this project; `D_{a,b} = Y^{(a,b)}_{(a,b)} = t^b − t^{b−1} + 1` is
paper-side (Theorem C at `m = b`), and in Lean that polynomial is simply written down. What this
session adds is exponent generality in the obstruction **mechanism**. Registry parent
`thm-D-product-form-obstruction` stays `proved`, not `lean-verified`.
