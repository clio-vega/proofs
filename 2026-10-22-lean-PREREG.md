# PRE-REGISTRATION — LEAN 2026-10-22, target: eps-profile determines the multiplicity set

Written BEFORE any `lake build` of the target module. Timestamp below is the shell's.

## Target (informal)

For finite `M, M' ⊆ ℕ₊` and `f_M(d) := #{a ∈ M : d ∣ a}`:
  (∀ d ≥ 1, f_M(d) = f_M'(d))  →  M = M'.

## Route I intend to take, declared in advance

NOT Möbius inversion. **Strong downward induction on the largest element of the
symmetric difference**: let `n = max (M ∆ M')`. Every multiple of `n` other than `n`
itself and `0` is `> n`, hence lies in `M ↔ M'` by maximality; so the fibre over `d = n`
differs by exactly the membership of `n`. Contradiction.

I predict Möbius resists because Mathlib's `ArithmeticFunction.moebius` inversion
lemmas are stated for sums over `Nat.divisors n` of a *fixed* `n`, and my sum is over
an arbitrary `Finset ℕ` support — the index sets do not line up without building a
bridge I do not need.

## Lesion register — PREDICTIONS

### L1. `M, M'` finite
**Prediction: IRREMOVABLE, and NOT A BINDER.** Finiteness lives in the *type*
(`Finset ℕ`), so there is no hypothesis to delete — grade it **TYPE-LEVEL / not
lesionable** unless I restate over `Set ℕ` with `Set.ncard`. If I do restate that way I
predict the theorem becomes **FALSE**, not merely unprovable: `Set.ncard` of an infinite
set is the junk value `0`, so `M = ℕ₊` and `M' = {even numbers}` both give every fibre
infinite, hence `ncard = 0` in every fibre, hence equal profiles with `M ≠ M'`. That
would be a junk-value *counterexample*, which is a different animal from a junk-value
*false removability* (10-21's `w≠0`).

### L2. `0 ∉ M` (i.e. `M ⊆ ℕ₊` rather than `M ⊆ ℕ`)
**Prediction: REMOVABLE.** This contradicts the brief, which says "expect this one to
bite". My reason: the brief is right that `d ∣ 0` for every `d`, so `0` sits in *every*
fibre — but that is exactly what exposes it. Take `d` larger than every nonzero element
of `M ∪ M'`; then the fibre over `d` is `(M ∪ M') ∩ {0}`, so the profile *reads off*
`0 ∈ M` directly. `0` is the one element the unbounded quantifier can isolate.
So I predict the theorem is true over all of `Finset ℕ` and the hypothesis deletes.

### L3. Quantifier range `∀ d ≥ 1`
**Prediction: SPLITS — carries a hypothesis in one restriction and not the other.**
- restricted to `d ∈ [1, max (M ∪ M')]`: **REMOVABLE / equivalent.** The witness `d = n`
  the proof uses satisfies `n ≤ max (M ∪ M')`.
- restricted to `d ∈ [1, max M]` (one side only): **IRREMOVABLE.** Counterexample
  `M = ∅`, `M' = {5}`: `max ∅ = 0`, range `[1,0] = ∅`, hypothesis vacuous, conclusion false.
- `d ≥ 1` vs `d ≥ 0`: `d = 0` adds the fibre `{a : 0 ∣ a} = {0}`, which is *information*,
  so allowing `d = 0` is a strictly weaker hypothesis to satisfy... no: allowing `d = 0`
  gives the prover MORE hypotheses, so it is a **weakening of the theorem**. Predict the
  `1 ≤ d` side condition is **REMOVABLE** from the hypothesis (I never use `d = 0`).

### L3'. Interaction prediction (the one I most want to test)
L2 and L3 are **not independent**. I predict that over `Finset ℕ` with `0` allowed, the
range restriction to `[1, max (M ∪ M')]` flips from REMOVABLE to **IRREMOVABLE**:
`M = {0}`, `M' = ∅` has `max (M ∪ M') = 0`, empty range, vacuous hypothesis, `M ≠ M'`.
So the two lesions are individually survivable and **jointly fatal**. If that is right,
testing them one at a time — which is what a lesion pass normally does — would grade
both REMOVABLE and miss the pair.

## Axiom check plan
`Lean.collectAxioms` on the main declaration, **plus a planted `sorry` positive control**
in the same module, because a sorry detector gave a false all-clear on 10-21 c1.

## Declaration census plan
Environment-derived (`#eval` over the environment, filtered by namespace), compared as a
SET in both directions against the source. Totals agreeing proves nothing.
Thu Oct 22 07:22:23 UTC 2026
