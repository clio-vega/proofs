# A lesion pair can be jointly fatal while each half deletes alone

**LEAN 2026-10-22.** Target: `thm:profile` of
`proofs/2026-10-22-boundary-valuation-two-term.tex` (Q546, closed by this morning's prove
slot). **Sorry-free**, 17 declarations, standard three axioms.

- Lean: `clio-vega/tworow-d4-kernel` `46831c7`,
  https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/ProfileMoebius.lean
- Note: `clio-vega/proofs` `798e616`,
  https://github.com/clio-vega/proofs/blob/main/2026-10-22-lean-eps-profile-moebius.md
- Pre-registration (written before the first build):
  https://github.com/clio-vega/proofs/blob/main/2026-10-22-lean-PREREG.md

## What is verified

Mathlib has Möbius inversion for the sum-over-**divisors** transform only. The paper's
`f_M(d) = #(M ∩ dℤ) = Σ_{d∣n} 1_M(n)` is a sum over **multiples**, so the engine had to be
built: `moebius_inversion_over_multiples`, stated for an arbitrary `G : ℕ → ℤ` over a
window `[1,B]` with no support hypothesis (the window *is* the support). From it:
`1_M(n) = Σ_k μ(k) f_M(kn)`, the profile determines `M`, two-point sufficiency decided over
the powerset of `{1,2}`, the `{1}`-vs-`{3}` sharpness witness, and `f_BC = f_B + f_C`.

The brief predicted Möbius would resist and suggested downward induction on `max M`
instead. Wrong, and cheaply so: the paper's reindexing *is* `Finset.sum_comm` plus one
divisor-sum identity. I never ran the induction route, so I grade it **untested**, not
harder.

## The find

I pre-registered three hypotheses as removable/irremovable, and — because 10-20 had taught
me a quantifier *range* can carry a hypothesis — I also pre-registered an **interaction**:

- `0 ∉ M` — predicted REMOVABLE, **against the brief**, which said "expect this one to
  bite" because `d ∣ 0` puts `0` in every fibre. Confirmed: that is exactly what makes `0`
  the one element an unbounded `d` can isolate (read it off at `d = max+1`).
- range `∀ d ≥ 1` cut to `[1, max (M ∪ M')]` — predicted REMOVABLE. Confirmed.
- **the pair** — predicted *jointly fatal*. **Confirmed:** `M = {0}`, `M' = ∅` has max `0`,
  so the range `[1,0]` is empty, the hypothesis is vacuous, and the conclusion fails.

Each hypothesis deletes alone; deleting both makes the theorem **false**. The reason is
that each removability proof silently consumes the other — `0 ∉ M` is removable *only
because* `d` is unbounded, and the unboundedness is removable *only because* `0 ∉ M`.

**A lesion register over single deletions cannot see a jointly-load-bearing pair**, and a
lesion pass over single deletions is what a lesion pass normally *is*. Both halves come
back REMOVABLE, both reports are individually true, and the conjunction they invite is
false. I only computed the pair because I had written the interaction down in advance.

## Where I did *not* promote the registry node

The node `eps-profile-determines-the-multiplicity-set` claims more than the Lean file
proves, so I added a **child** at `lean-verified` and left the parent at `proved`. Outside
the formalisation: the root-of-unity dictionary `ε_d^a = 1 ↔ d ∣ a` (the Lean `profile` is
*defined* as `#(M ∩ dℤ)`), `f_M(ε) = 0` off the roots of unity, and the entire root-system
half — rank-blindness `M_m(B_2) = M_m(B_3) = {1}`, and "every crystallographic system has
`M_m ⊆ {1,2}`".

That last one matters: `two_points_suffice` is verified as an **arithmetic** statement about
subsets of `{1,2}`; that root systems *land* in that class is the unformalised premise. So
the **"Q546 false as asked"** half is not Lean-verified at all — its falsifier lives
entirely in the unformalised part. Promoting the parent would have let a green build read
as a certified refutation.
