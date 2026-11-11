# LEAN 2026-11-11: Q565 `thm:red` (a reciprocal pair is a restriction, not a vanishing)

**Project:** `lean/tworow_d4_kernel`, new module `TworowD4Kernel/ReciprocalRestriction.lean`.
Commit `bbdbb83` was pushed to `origin/main`. I checked it with `git show --stat` (1 file, +461) and
`git status -sb` reads `## main...origin/main`.
**Paper:** `proofs/2026-10-24-Q565-zero-locus-in-M.tex`, `thm:red`, L241–L291.

## What builds sorry-free

`lake build` returned rc=0. I read the status directly, not through a pipe. The module was built
(`✔ [3481/3482]`) with no warnings and 0 `sorry`.

| Layer | Declarations | Status |
|---|---|---|
| Cross-factor identity (primary target, no `thm:alt`) | `crossFactor_mul_crossFactor_of_recip` (b), `crossFactor_eq_one_of_self_recip` (a) | **sorry-free, unconditional** |
| `B(T)` layer | `Bfactor_eq_zero_of_meets_pair`, `Bfactor_erase_pair`, `Bfactor_eq_zero_of_mem_self_recip`, `Bfactor_erase_self` | **sorry-free, unconditional** |
| Alternant layer (`psiAlt`) | `psiAlt_erase_pair`, `psiAlt_erase_self` | **sorry-free, unconditional** |
| `eq:Psidef` (`psiHead`) | `thm_red_b`, `thm_red_a` | **sorry-free, CONDITIONAL** on `AltExpansion` (`thm:alt` as a `Prop` binder) at exactly two sets: `M` and the reduced set |

The cross-factor identity holds under exactly three hypotheses: `w_i ≠ 0`, `w_i ≠ w_a` and
`w_i ≠ w_b`. It does not need `w_a ≠ 0` as a hypothesis, because that follows from `w_a w_b = 1`.

There are no sorries. The one thing left unformalised is `thm:alt` itself. It is still the gap for
the whole development, and it enters only as a hypothesis. The control `altExpansion_rank_zero` shows
the hypothesis can be satisfied, but only at `r = 0`. That is a non-vacuity check and is **not**
evidence for `thm:alt` at `r ≥ 1`.

## `#print axioms`

All 15 headline declarations, including `thm_red_a`, `thm_red_b`, `psiAlt_erase_pair` and the
cross-factor identities, give `[propext, Classical.choice, Quot.sound]`. In the same invocation, the
planted `theorem bad_control : (1:ℕ) = 2 := by sorry` gave `[sorryAx]`.

## Carrier

I used a `Finset ι` of **indices** with weights `w : ι → K`, which is the paper's `conv:psi`. The
removal `M ∖ {w_a, w_b}` is `(M.erase a).erase b`, so exactly one index is removed for each element
of the pair.

## Lesions

| # | Prediction | Verdict |
|---|---|---|
| L1 drop `a ≠ b` | brief: refutable; mine (pre-build): collapses to (a) | **REMOVABLE.** The brief's prediction failed. Lean *proves* the lesioned (b) (`L1_psiAlt_erase_pair_without_ne`). |
| L2 `w_a w_b=1 ↦ w_a=w_b` | vacuous under admissibility, refutable without | **CONFIRMED.** The witness is `w=(2,2,3)` with cross factors `5/2·5/2 ≠ 1`. |
| L3 `Finset` of values | brief: refutable (repeated entry) | **ABSENT.** `L3_value_carrier` is a corollary. |
| {L1,L3} | brief: refutable | **ABSENT.** `L1_L3_value_carrier_without_ne` proves it. |
| L4 drop admissibility | refutable at the `B` layer | **CONFIRMED.** Witness `w=(2,1/2,2)`, where the repeated value `w_2=w_0` gives `x/0=0`, so `0 ≠ 3/4`. |

**Finding:** the Q565-A carrier lesson does not carry over. In `prop:recip` the collapse
`{a,a}={a}` lost a square inside a **product**. In `thm:red` the set `{a,b}` is only ever
**removed**, and when `a = b` removing one copy is the correct operation. The repeated entry the
brief expected is real (`wL4`), but admissibility excludes it (L4), not the choice of carrier.

## Registry

In `rick-q488-converse.json` I added a new child, `clio-Q565-reduction-alternant-layer`, with trust
`lean-verified` and the six unconditional declarations. Its parent, `clio-Q565-reduction-theorem`,
**stays `proved`**, because its Ψ statement rests on `thm:alt`. This follows the precedent of
`clio-Q559-recip-pair-structural-vanishing-B-layer`. `registry_validate` gives OK. As a negative
control, a fault I planted (the missing note file) was flagged.

No correction to the `.tex` is needed.
