# Lean: the migration potential — `Δφ = 1` grades, `Δφ = ±1` does not

**Date:** 2026-10-06 (LEAN c1)
**Project:** `projects/lean/tworow_d4_kernel`
**File:** `TworowD4Kernel/MigrationPotential.lean` (namespace `TworowD4Kernel.MigrationPotential`)
**Registry:** `proofs/registry/migration-length-grading.json`
**Paper proofs formalised:** `proofs/2026-10-05-c3-migration-length-grading.tex` (Theorem 1),
`proofs/2026-10-06-migration-vertical-step-gap.tex` (Corollary 2, `cor:irred`)
**Source:** Purbhoo, *Puzzles, tableaux and mosaics*, arXiv:0705.1184 — §3.1 (step rule), §4 (the wake)

## Status: sorry-free, 6 declarations, 0 sorries

`lake build` exit 0 (`PIPESTATUS[0]=0`), both for the single module and for the whole
project after adding the import to the root aggregator (3193 jobs, 189 s).

| Declaration | Content | Axioms |
|---|---|---|
| `length_eq_potential_diff` | Task 1 — telescoping: `Δφ = 1` on every step ⇒ `l.length = φ b − φ a` | `propext, Quot.sound` |
| `length_eq_of_endpoints` | Task 2(a) — gradedness: two chains, same endpoints ⇒ same length | `propext, Quot.sound` |
| `sum_pow_length` | Task 2(b) — `∑_{M∈F} q^{ℓ M} = ‖F‖·q^{ℓ₀}` when `ℓ` is constant on `F` | `propext, Quot.sound` |
| `no_potential_of_step_symm` | Task 3 / `cor:irred` abstract half — a mutually inverse pair admits no `+1` potential | `propext, Quot.sound` |
| `pm_one_not_graded` | Task 3 — the `±1` relaxation: explicit witness, chains of length 1 and 3, same endpoints | `propext, Quot.sound` |
| `pm_one_witness_has_no_potential` | the two above joined: the witness relation admits **no** `+1` potential, and for the `cor:irred` reason | `propext, Quot.sound` |

No `Classical.choice` — a **strict subset** of the permitted three.

`sorry` scan: recursive (`grep -rn --include=*.lean -w sorry . --exclude-dir=.lake`),
**validated against a planted canary** — 11 hits before, 12 with the canary appended,
11 after restore (`cmp` byte-identical). All 11 standing hits are prose in docstrings
("`sorry`-free", "No `sorry`"), none a tactic. The authoritative check is the axiom table
above: no `sorryAx` anywhere.

## What Mathlib already had

**Nothing matching.** This was the mandatory pre-flight, and I searched the *shape*:

* `Mathlib/Data/List/Chain.lean` — no conclusion of the form
  `l.length = f (getLast …) − f (head …)`. There is `IsChain.backwards_induction` with a
  `getLast` hypothesis (the right induction skeleton) but no arithmetic conclusion.
  `grep` for `telescop` across all of Mathlib: no chain-length result.
* `Mathlib/Order/Grade.lean` — `GradeOrder` is the nearest object, and it is **not** this
  statement. It fixes the relation to the covering relation `⋖` of a preorder and the grade
  to land in a graded order, then *derives* `grade b = grade a + 1` from `a ⋖ b`. Here the
  implication runs the other way: the `+1` is the **hypothesis** and `step` is arbitrary.
  A migration step is not presented as a covering relation of any order on mosaics, so
  `GradeOrder` does not apply even after transport.

So the telescoping is mine. It is also three lines, which is presumably why it is not in
Mathlib — and that is the honest size of Task 1.

**Toolchain note for the next slot:** in this Mathlib (`v4.30.0`) `List.Chain r a l` is
gone; the object is Batteries' `List.IsChain R l` on a *single* list, with
`List.isChain_cons_cons : IsChain R (a :: b :: l) ↔ R a b ∧ IsChain R (b :: l)` as the
destructor and `IsChain.singleton`/`cons_cons` as constructors. The brief's statement
shape (`List.Chain step a l`) does not typecheck.

## What the `±1` counterexample shows, and what it does not

The 10-05 note's Theorem 1 rests on `Δφ = 1` at **every** migration step. The gap
(registry `gap-vertical-step-orientation`): the strictly vertical 180° step exists in two
directions, `Δφ = ±1`; **91 of 707** computed steps were strictly vertical and all 91 took
`+1`; today's PROVE slot widened the null to **3111** steps with a `−1` candidate
*available* in 0 of them. There is still no proof.

`pm_one_not_graded` converts that silence into a type-checked statement of what is at
stake: weaken the hypothesis to `±1` and the conclusion **genuinely fails** — the witness
is `step x y ↔ y = x ± 1` on `ℤ`, `φ = id`, and the chains `0 → 1` (length 1) and
`0 → 1 → 2 → 1` (length 3), same endpoints, different lengths. So the `+1` is
**load-bearing, not cosmetic**.

`pm_one_witness_has_no_potential` makes the pair a real dissociation rather than an
accident: the witness relation admits no `+1` potential at all, for exactly the reason
`cor:irred` names — `0 → 1` and `1 → 0` are both steps. The counterexample is not an
artefact of a bad choice of `φ`.

**It does not show the gap is real.** Nothing here bears on Hypothesis (G). The Lean file
is entirely abstract; it says *if* the backward vertical direction is realised, the
`ℓ`-grading dies and no better height function saves it. Whether it is realised is a
question about mosaics and is untouched.

## Honest scope line — what is NOT formalised

Everything geometric. In particular:

1. **That `φ(a,b,c,e) = b + e` has `Δφ = 1` on actual migration steps.** This is the
   exhaustive classification of the 10-05 note (90 hexagon regions × every tiling ×
   every rhombus × every rotation, under Purbhoo's constraints (C1)–(C3)). Registry node
   `classification-of-steps`, trust `computed`, and it stays `computed`.
2. **Proposition 1 of the 10-06 note** — that the two strictly vertical steps are the two
   tilings of the zonogon `R ⊕ [0,1]u₆₀`, interchanged by its central symmetry, hence
   mutually inverse. Another exhaustive enumeration. It is the **hypothesis** of
   `no_potential_of_step_symm`, supplied from outside.
3. **Hypothesis (G)** (`hyp-G-no-strip-below`) and **(G2)** (`hyp-G2-one-rhombus-per-hexagon`).
   Both `in-progress`, both unproved, both untouched by this file.
4. **The closed form for `ℓ₀`** (`closed-form-ell0`) — checked on 361 fibres and 487
   journeys, `computed`, not formalised.

So: Theorem 1 of the 10-05 note is **not** Lean-verified, and the parent registry node is
**not** promoted. `unproved` is not `unformalised`, and the reverse holds too — what is
formalised here is the abstract skeleton, and the skeleton was never the doubtful part.

## Two slots, one argument — NOT independent corroboration

The brief required that this slot and today's PROVE slot not be allowed to agree by
construction. They **do** agree, and the agreement carries less information than it looks
like:

* PROVE's Corollary 2 proof: the two vertical steps are mutually inverse, so
  `Δψ(T₊→T₋) = −Δψ(T₋→T₊)` for any `ψ`, and both cannot be `+1`.
* `no_potential_of_step_symm`: `step a b`, `step b a`, `ψ b = ψ a + 1`, `ψ a = ψ b + 1`,
  `omega`.

**That is the same argument, type-checked.** Different mechanism for the *geometry*
(PROVE: exhaustive enumeration over hexagon regions; Lean: none — it is a hypothesis), but
identical mechanism for the *implication*. Logging it this way rather than banking two
confirmations: `a-second-instrument-can-be-the-first-instrument`.

What the Lean version does add is a **separation of the two halves**. PROVE's corollary
bundles a geometric claim with a one-line algebraic deduction; Lean forces the geometric
claim out into a hypothesis and shows the deduction needs nothing else — no finiteness, no
connectivity, no properties of mosaics. That tells you exactly where the remaining work is,
which is the one thing a prose corollary of that shape tends to blur.
