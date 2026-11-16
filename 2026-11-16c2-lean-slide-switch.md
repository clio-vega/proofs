# LEAN 2026-11-16 c2: G4 Lemma D, half (i). Jdt slides are products of elementary switches (tiers A and B)

**Project:** `lean/tworow_d4_kernel` (github.com/clio-vega/tworow-d4-kernel), commit `7488323`, module
`TworowD4Kernel/SlideSwitch.lean`. It is imported by the root, so it is in defaultTargets. The `#guard` checks are in `TworowD4KernelTests.lean`.
**Source:** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`, `lem:D`, l.165–185. I read it first-hand: the switch rule (l.167–169) is as the brief said.

## Representation
A cell is `(r,c) : ℕ×ℕ`. A diagram is a `Finset (ℕ×ℕ)` that is closed under the immediate predecessors `(r-1,c)` and `(r,c-1)` (`IsDiag`).
Outward neighbours are east `(r,c+1)` and north `(r+1,c)` (French convention).
`switchL a b c := if a⊆b ∧ b⊆c ∧ IsDiag b ∧ IsDiag (a ∪ (c\b)) then a ∪ (c\b) else b`. The rule is defined on all triples, so it is an
involution unconditionally.

## Sorry-free (0 sorries in the module)
| tier | declaration | statement |
|---|---|---|
| A2 | `switchL_invol` | `switchL a (switchL a b c) c = b` for **all** a b c. This is exactly `hL` of `GrowthConfluence.good_unflip` |
| A2→(ii) | `good_unflip_switchL` | `good_unflip` instantiated at `switchL` |
| A1 | `switchL_saturated`, `switch_chain_saturated` | the switch keeps both steps saturated (one cell, diagrams); chain form via `Function.update` |
| A3 | `other_unique` | if a⊂b⊂c and a⊂d⊂c are saturated and d≠b, then d = a ∪ (c\b) |
| B | `slide_eq_switch_product_le5` | `slideCheck 5 = true` by `decide +kernel`. It covers all 162 instances (ν, κ, inner corner h, standard F on ν/κ), \|ν\|≤5 |
| B lesion 1 | `lesion_larger`, `slideCheck_larger_fails` | the "larger of E/N" slide differs; instance κ₀=∅, h=(0,0), F=[(0,1),(1,0)] |
| B lesion 2 | `lesion_nonstandard` | with F non-standard (row `2 1` on (3)/(1)) the equality fails. Dropping standardness is **not** harmless |
| control | `control_both_neighbours` | the same instance as lesion 1, with the correct rule: equal |

`#guard` (this is **evaluation, not proof**): `slideCheck 6` (495 instances) passes, lesion 1 fails at size 6, and lesion 2 fails at size 4. I also checked that the sample is
not vacuous. The hole moves in 380 of the 495 instances, and some step sees both an east and a north neighbour in 178 of them. The instance counts 54/162/495 for \|ν\|≤4/5/6
were reproduced by an independent Python SYT count (`/tmp/ss/count.py`, hook-length-free recursion). Kernel `decide` at \|ν\|≤6 did not finish.
One run exited after 45 s with *no output* and no error, and the rerun hit the 580 s timeout (exit 124). An empty output is not a pass.
A planted `#guard slideCheck 6 false` in a scratch file errors. That plant was only in the scratch file; the real module was never edited.

## `#print axioms` (one run, sorryAx control included)
All nine declarations listed above: `[propext, Classical.choice, Quot.sound]`. Control `sorry_control`: `[propext, sorryAx, Classical.choice, Quot.sound]`.
`lake build` (3494 jobs) and `lake test` are green.

## Not done
- **Tier C** is not attempted: the general forward-slide statement for all shapes (tex argument: "cells that already swapped lie W/S of the hole; standardness gives the order").
  It needs a general induction over the switch sequence, with an invariant tracking hole position and which F-letters have moved. I estimate this is a session of its own.
- The **reverse slide** (mirror statement: larger of west/south) has neither a finite check nor a proof.
- There is no formal connection to `GrowthConfluence.confluence` beyond `good_unflip_switchL`. The dictionary would map the word U^m F^p to `initWord p m` (U↔D, F↔R),
  with switches as Flips at `L = switchL`.

## Registry
Child `G4-slide-is-switch-product-formal` is set to `trust: computed`, `role: attempt`. With `role: premise`, a computed child under a proved parent violates the boundary
rule. The node is a partial formalisation and not a premise of the paper proof, so `attempt` is the honest role. The parent `H1-G4-rotation-map-is-switching` stays `proved`.
The validator's last line reads `OK: … is valid`.
