/-
Copyright (c) 2026 Clio. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Clio
-/
import Mathlib.Algebra.Field.ZMod
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Algebra.Polynomial.SpecificDegree
import Mathlib.Tactic.ComputeDegree

/-!
# The 𝔽₂-irreducibility of `X² + X + 1`

This file formalises the load-bearing finite-field input to the **unique 2-adic root**
result for the two-row `d = 4` law in the residue classes `b ≡ 2, 3 (mod 4)`
(`~/projects/proofs/2026-06-09-tworow-d4-b23-2adic-single-candidate`).

That result reduces, modulo `2`, to
`Q_b ≡ m · (m² + m + 1)^{⌊(b-1)/4⌋} (mod 2)`,
where `m` (here the indeterminate `X`) is a **simple** root of `X³ + X² + X` over `𝔽₂` and
`X² + X + 1` is **irreducible** over `𝔽₂`. These are exactly the two hypotheses Hensel's lemma
needs to conclude that `Q_b` has at most one root in `ℤ₂`, hence at most one rational root.

## Main results

* `TworowD4Kernel.irreducible_X_sq_add_X_add_one`: `X² + X + 1` is irreducible over `ZMod 2`.
  Proof: it has degree `2` and no root in the two-element field, so the degree-`≤ 3`
  no-root criterion `Polynomial.irreducible_of_degree_le_three_of_not_isRoot` applies.

* `TworowD4Kernel.not_X_dvd_X_sq_add_X_add_one`: `X` does not divide `X² + X + 1` over
  `ZMod 2`. Equivalently `X` is a **simple** root of `X · (X² + X + 1) = X³ + X² + X`: it has
  nonzero constant term, so `X ∤ (X² + X + 1)`. This is the simple-root input to Hensel's lemma.

## References

The informal proof is `2026-06-09-tworow-d4-b23-2adic-single-candidate`; the irreducibility and
simple-root facts are the finite-field hypotheses of the Hensel uniqueness of `ρ_b`.
-/

namespace TworowD4Kernel

open Polynomial

/-- The natural degree of `X² + X + 1` over `ZMod 2` is `2`. -/
theorem natDegree_X_sq_add_X_add_one : (X ^ 2 + X + 1 : (ZMod 2)[X]).natDegree = 2 := by
  compute_degree!

/-- Neither element of `ZMod 2` is a root of `X² + X + 1`: evaluating gives `0² + 0 + 1 = 1`
and `1² + 1 + 1 = 1`, both nonzero. -/
theorem not_isRoot_X_sq_add_X_add_one (x : ZMod 2) :
    ¬ (X ^ 2 + X + 1 : (ZMod 2)[X]).IsRoot x := by
  simp only [IsRoot.def, eval_add, eval_pow, eval_X, eval_one]
  fin_cases x <;> decide

/-- **Irreducibility input to Hensel uniqueness.**

`X² + X + 1` is irreducible over the two-element field `ZMod 2`. Since it has degree `2` and no
root in `ZMod 2` (neither `0` nor `1` is a root), it is irreducible by the degree-`≤ 3` no-root
criterion. -/
theorem irreducible_X_sq_add_X_add_one :
    Irreducible (X ^ 2 + X + 1 : (ZMod 2)[X]) := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  apply Polynomial.irreducible_of_degree_le_three_of_not_isRoot
  · rw [Finset.mem_Icc, natDegree_X_sq_add_X_add_one]
    exact ⟨by norm_num, by norm_num⟩
  · exact not_isRoot_X_sq_add_X_add_one

/-- **Simple-root input to Hensel uniqueness.**

`X` does not divide `X² + X + 1` over `ZMod 2`. Since `X · (X² + X + 1) = X³ + X² + X`, this says
exactly that `X` is a **simple** root of `X³ + X² + X` (multiplicity `1`): the cofactor
`X² + X + 1` does not vanish at the root, which here is visible as its nonzero constant term
`(X² + X + 1).coeff 0 = 1 ≠ 0`. -/
theorem not_X_dvd_X_sq_add_X_add_one : ¬ (X ∣ (X ^ 2 + X + 1 : (ZMod 2)[X])) := by
  rw [X_dvd_iff]
  simp only [coeff_add, coeff_X_pow, coeff_X, coeff_one_zero]
  decide

end TworowD4Kernel
