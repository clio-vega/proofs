# LEAN 2026-10-23 — the ρ-shifted vanishing law: index arithmetic + equal-column determinant

**Session:** Lean (one target). **Status: sorry-free.**

## Target

Forward direction of `thm:vanish`,
`proofs/2026-10-22-c2-Q553-residue-groupings-and-head-collapse.tex` (L388–L397).

* **Project:** `projects/lean/tworow_d4_kernel` (Mathlib `v4.30.0`, Lean `v4.30.0`)
* **New module:** `TworowD4Kernel/RhoShiftVanishing.lean` (214 lines, 18 declarations)
* **Main declaration:**
  `TworowD4Kernel.RhoShiftVanishing.det_alternant_eq_zero_of_rhoTie_lt`

```lean
theorem det_alternant_eq_zero_of_rhoTie_lt {K : Type*} [Field K] (x : Fin r → K)
    (θ : Fin r → ℤ)
    (h : ∃ j k : Fin r, j < k ∧ θ j - θ k = (k.val : ℤ) - (j.val : ℤ)) :
    (Matrix.of fun i j => x i ^ (rhoShift θ j)).det = 0
```

The module is reached by the `TworowD4Kernel.+` glob in `lakefile.toml`; the root
aggregator was not edited.

## What builds, and what it does not cover

Everything in the module builds **sorry-free**. `grep -n "sorry\|admit\|axiom "` over the
source returns nothing. Full-project `lake build`: 3479 jobs, completed; every one of the
81 source modules has an olean (set comparison, both directions — 5 oleans exist with no
source, stale scratch modules predating this session).

The brief asked for three pieces, deliberately not grouped. All three are present:

| | declaration | content |
|---|---|---|
| **A** | `exists_rhoTie_iff_not_injective` | `(∃ j ≠ k, θ j − θ k = k − j) ↔ ¬Injective (rhoShift θ)` |
| **B** | `det_zpow_eq_zero_of_not_injective`, `det_npow_eq_zero_of_not_injective` | two equal columns ⟹ `det = 0` |
| **C** | `det_alternant_eq_zero_of_rhoTie`, `det_alternant_eq_zero_of_rhoTie_lt` | the composition; `_lt` is the paper's hypothesis verbatim |

**(A) is an `↔`, not merely sufficiency.** The paper's proof says "is exactly"; Lean
holds it to that, and the reverse direction is what the converse (Q559) consumes.

**The layer this does NOT cover — and it is the honest headline.** The paper's conclusion
is `Ψ^C_{(1^r),θ}(M) = 0`. What is formalised is `det(x_i^{m_j}) = 0`. The passage from
one to the other is the alternant expansion `thm:alt` (L341–L363), which is **not
formalised**. A Lean-verified summand is not a Lean-verified sum. Accordingly the registry
node `clio-Q553-rho-shift-vanishing-law-forward` **stays at `proved`**; only the new child
`clio-Q553-rho-shift-vanishing-determinant-layer` is promoted to `lean-verified`.

## Where I expected the content, and where it was

The brief predicted the content would be in the index bridge (A) — "the `−1` is inert,
the `Fin r → ℤ` coercion is not" — and pre-registered that expectation. **The prediction
was wrong in a specific, cheap way.** `rhoShift_eq_iff` is `unfold; omega`, one line, and
the coercion cost nothing: `omega` handles `(j.val : ℤ)` natively. The reason the
indexing convention could not bite is structural and worth recording:

> The paper is 1-indexed and `Fin r` is 0-indexed, so the paper's `m_j = θ_j + j − 1` is
> this file's `rhoShift` **plus 1**. But the *hypothesis* only ever mentions the
> **difference** `k − j`, which is convention-independent
> (`(k+1) − (j+1) = k − j`); and the *conclusion* is invariant under a constant shift of
> all exponents. So the convention cannot reach either end of the statement.

That is `injective_rhoShiftPaper_iff` (via `injective_add_const_iff`), which makes the
inertness a Lean-checked fact rather than a remark — this is L4, below.

## Pre-registered lesions, graded

| | lesion | prediction | **Lean grade** |
|---|---|---|---|
| **L1** | `j ≠ k` weakened from `j < k` | REMOVABLE / equivalent | **REMOVABLE — confirmed.** `rhoTie_lt_iff_ne` proves the two hypotheses *equivalent* (`lt_or_gt_of_ne` + `omega` for the swapped equation), and `det_alternant_eq_zero_of_rhoTie_lt` is derived from the `≠` form by `ne_of_lt`. |
| **L2** | `x i ≠ 0` | SPLIT: irremovable for `zpow`, absent for `npow` | **ABSENT in BOTH arms — prediction wrong. This is the session's find.** See below. |
| **L3** | `r ≥ 2` / `2 ≤ Fintype.card` | ABSENT (the quantifier range carries it) | **ABSENT — confirmed, and located.** `rhoTie_forces_two_le` *derives* `2 ≤ r` from the hypothesis `∃ j k : Fin r, j ≠ k`. The binder is absent, not relocated: there is no upstream lemma stating it, because the existential's own range implies it. |
| **L4** | the `−1` in `m_j = θ_j + j − 1` | REMOVABLE (any constant shift works) | **REMOVABLE — confirmed, and strengthened.** `injective_add_const_iff` proves it for an arbitrary constant `c` over an arbitrary domain, so the whole family of conventions is covered at once, not just `−1`. |

### L2 is the find: invertibility is not where I filed it

I predicted `x i ≠ 0` would be **irremovable for the `zpow` statement** — the Laurent
alternant is the one the paper actually claims (Setup L41 specialises the `w`-values to
nonzero scalars), and negative exponents looked like they must demand units.

Lean says **ABSENT**. `det_zpow_eq_zero_of_not_injective` is proved over a bare `Field`
with no hypothesis on `x` at all, by the *same four-line proof* as the `npow` arm. Two
reasons, and both were in front of me:

1. Mathlib's `zpow` is **total** on a `GroupWithZero` (`0 ^ (-1) = 0`), so the matrix
   entries exist whether or not `x i` is a unit;
2. `Matrix.det_zero_of_column_eq` needs two columns to be **equal**, and nothing more.
   Equality of columns is insensitive to whether the entries are invertible.

So the Lean statement is strictly *stronger* than the paper's (it covers `x i = 0`), and
still true — a faithful generalisation, not a weakening.

**Where invertibility actually lives.** It is spent in exactly one place, and I proved
that place so the hypothesis would have somewhere to be rather than merely being absent:

```lean
theorem det_zpow_shift {K : Type*} [Field K] (x : Fin r → K) (hx : ∀ i, x i ≠ 0)
    (m : Fin r → ℤ) (c : ℤ) :
    (Matrix.of fun i j => x i ^ (m j + c)).det
      = (∏ i, x i ^ c) * (Matrix.of fun i j => x i ^ m j).det
```

This is the monomial-rescaling bridge — the brief's own route from the `zpow` statement
down to a shifted `npow` one — and it is where `zpow_add₀` demands `x i ≠ 0`. **The
vanishing law does not use it.** I had reasoned "the `zpow` statement needs units"
because the *bridge to it* needs units; the statement itself does not. The hypothesis
belongs to the theorem that rescales, not to the theorem that vanishes — the same shape
as the 10-19 note that a genericity hypothesis belongs to the theorem that separates.

**Which statement the `.tex` was claiming:** the `zpow` one. `m_j ∈ ℤ` and the `w`-values
are nonzero scalars, so the paper's alternant is Laurent. The `npow` version is proved
*separately* (`det_npow_eq_zero_of_not_injective`) and is not substituted for it.

### Free lesion pass — the `unusedVariables` linter

Read, as instructed, and **it flagged nothing in this module** (the only two warnings in
the 3479-job build are pre-existing, in `ParityObstruction` and `IotaAdjunction`). No
binder was silently non-load-bearing. I did not underscore anything to quiet the linter.

## Controls

**Planted-sorry control — FIRED.** A planted `theorem planted_sorry_control : 2 + 2 = 5 := by sorry`
produced `warning: TworowD4Kernel/RhoShiftVanishing.lean:173:8: declaration uses 'sorry'`,
then was removed. Worth recording precisely: the build *still printed*
`Build completed successfully`. A `sorry` is a warning, not an error — so the build's
success line carries **no information** about sorry-freeness, in either direction. The
grep and the warning line are the measurement.

**Positive control, vanishing arm.** `θ = (0, −1)` at `r = 2`, `x = (2, 3)` over `ℚ` —
the degenerate flipped term `Ψ^C_{(1,1),(0,−1)} ≡ 0` that was the discovery instrument on
10-22 (L453).
* `control_tie_hypothesis` — the tie hypothesis genuinely holds here (`decide`).
* `control_tie_det_zero` — and the determinant *is* `0`, computed directly through
  `Matrix.det_fin_two`, **not** by invoking the theorem. An independent reading.

**Positive control, non-vacuity arm.** `θ = (0, 0)`, same `x`.
* `control_no_tie_injective` — the ρ-shift is injective here, so the theorem does not
  apply (`decide`).
* `control_no_tie_det_ne_zero` — and the determinant is `1/2 − 1/3 = 1/6 ≠ 0`.

The second arm is the one that matters: it shows the vanishing law is not vacuously true
of every `θ`, which no amount of confirming the first arm could establish.

## Declaration census — set comparison, both directions

Environment scan (`run_cmd` over `env.constants`, non-internal, namespace-filtered)
against a source grep, compared **as sets**:

* source declarations: **18**
* environment declarations: **18** (19 minus the auto-generated `rhoShift.eq_1`)
* `comm -23` (source ∖ env): **empty**
* `comm -13` (env ∖ source): **empty**
* anonymous `example`s in source: **0** — so nothing is invisible to the environment scan

The count agreeing is not the evidence; the empty symmetric difference in both directions
is. A bare regex has over- and under-counted on different days, so the *sign* of a count
mismatch carries no information.

## `#print axioms`

All thirteen graded declarations depend on the standard three or fewer:

```
det_alternant_eq_zero_of_rhoTie_lt    [propext, Classical.choice, Quot.sound]
det_alternant_eq_zero_of_rhoTie       [propext, Classical.choice, Quot.sound]
exists_rhoTie_iff_not_injective       [propext, Classical.choice, Quot.sound]
det_zpow_eq_zero_of_not_injective     [propext, Classical.choice, Quot.sound]
det_npow_eq_zero_of_not_injective     [propext, Classical.choice, Quot.sound]
det_zpow_shift                        [propext, Classical.choice, Quot.sound]
injective_rhoShiftPaper_iff           [propext, Classical.choice, Quot.sound]
control_tie_det_zero                  [propext, Classical.choice, Quot.sound]
control_no_tie_det_ne_zero            [propext, Classical.choice, Quot.sound]
control_no_tie_injective              [propext, Classical.choice, Quot.sound]
rhoTie_forces_two_le                  [propext, Quot.sound]
rhoTie_lt_iff_ne                      [propext, Quot.sound]
control_tie_hypothesis                [propext]
```

No `sorryAx`, no custom axioms.

## Registry

* `clio-Q553-rho-shift-vanishing-law-forward` — **stays at `proved`**. Its conclusion is
  about `Ψ^C`; the alternant expansion is unformalised.
* **new child** `clio-Q553-rho-shift-vanishing-determinant-layer` — `trust: lean-verified`,
  `lean: RhoShiftVanishing.det_alternant_eq_zero_of_rhoTie_lt`.
* The unformalised layer already *is* its own node — `clio-Q553-alternant-expansion`
  (`proved`) — so no new node was needed for it; the dependency is named explicitly in
  both the parent's and the child's `approach`.

Registry file: `proofs/registry/rick-q488-converse.json`
(backup `…json.bak-1023-lean`).
