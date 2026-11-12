# LEAN 2026-11-12: Q559 `lem:reflect` is sorry-free

**Target.** `lem:reflect` from `proofs/2026-10-23-Q559-vanishing-converse.tex`, L304–L345:
ι Ψ^C_{(1^r),θ}(M) = (−1)^r e^{2r} Ψ^C_{(1^r),θ∨}(M), with m∨ = −m.
**Project.** `lean/tworow_d4_kernel` (github.com/clio-vega/tworow-d4-kernel), module
`TworowD4Kernel/Reflection.lean`, commit `c4cd1b5` (pushed, confirmed upstream with `gh api`).

## Statement as formalised
```lean
theorem lem_reflect (w : ι → K) (M : Finset ι) (r : ℕ) (θ : Fin r → ℤ)
    (hadm : Admissible w M) :
    psiHead (fun t => (w t)⁻¹) M r θ =
      (-1) ^ r * (∏ t ∈ M, w t) ^ (2 * r) * psiHead w M r (thetaVee θ)
-- thetaVee θ j = -θ j - 2 * j   (0-indexed j; the paper's -θ_j - 2j + 2 with 1-indexed j)
```
`psiHead` is the head-injection sum eq:Psidef. The identity holds at every point of any field
where `w` is nonzero and pairwise distinct on `M`, which specialises the identity in ℚ(w).

## Which θ∨: measured before writing the statement
`scratch/1112_lem_reflect.py` uses exact Fractions straight from eq:Psidef, at 60 random admissible
signed points with r ≤ 3 and n ∈ r…r+2:
- the paper's reading (θ∨_j = −θ_j − 2j, 0-indexed) holds at **60/60** points;
- the repository-`rhoShift` reading (+2) holds at only **10/60**. The 10 are exactly the points where the LHS is 0. This reading is
  refuted in Lean by `lesionVee_refutable` (r=1, M={0}, w=2, θ=0: −3 ≠ −12).

## Route: a deviation, recorded
The paper applies ι to `thm:alt`. The Lean proof applies the same three elementary identities
(L320–L325) to each head injection of eq:Psidef:
- each diagonal factor contributes −X_j^{2+2j} and shifts θ to θ∨;
- each pair j<k contributes X_j²;
- each cross factor contributes w_t².
Every X_j's exponents add up to 2r, and ∏X·∏_{t∉h} w_t = e. So **`lem_reflect` does not use `thm_alt`**. The computation is the
paper's, done on a finer index set. It is not new mathematics.

## Sorries
None in this module. `lake build` rc=0 (3484 jobs).

## Axioms (same invocation as a `sorryAx` control)
```
'bad_control_1112' depends on axioms: [sorryAx]
'TworowD4Kernel.Reflection.lem_reflect' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.Reflection.lem_reflect_of_ne_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.Reflection.lesionN_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.Reflection.lesionVee_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.Reflection.rhoShiftPaper_thetaVee' depends on axioms: [propext, Quot.sound]
```

## Lesions (each asked as "derive False", with predictions registered first)
- **N (drop w ≠ 0, which is equivalent to e ≠ 0): refutable.** `lesionN_refutable`, with r=1, M={0}, w=0,
  θ=0. The LHS is 1 (because 0⁻¹=0 and 0⁰=1) and the RHS is 0.
- **D (drop distinctness): removable in Lean.** `lem_reflect_of_ne_zero`. This is only the
  x/0=0 totalisation: when two values coincide, the denominator 1−a/b dies on both sides and both sides are junk 0.
  The statement is wider, not stronger, because coincidences do not occur in ℚ(w).
- **{N, D}: jointly fatal**, already through N's witness (|M|=1, so D is vacuous there).
- **θ∨ = repository reading: refutable** (see above).

## Registry
`clio-Q559-reflection-identity` in `registry/rick-q488-converse.json` went from `peer-reviewed` to `lean-verified`,
with `lean: Reflection.lem_reflect` and a scope note (proofs commit `2701246`). `registry_validate`:
OK. Its positive control on a corrupted copy reported 97 problems, but most of those are path artefacts of
running it from /tmp. So that control fired without being specific to this field.

## Not done
Stretch target Lemma 5 (order ≤ −r at w_n → ∞) was not attempted.
