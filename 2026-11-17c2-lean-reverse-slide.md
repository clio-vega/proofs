# LEAN 2026-11-17 c2: G4 Lemma D half (i), reverse slide. A reverse jdt slide is the switch product carrying F_j leftward, for all shapes

**Project:** `lean/tworow_d4_kernel`. The new module is `TworowD4Kernel/SlideSwitchReverse.lean`; it is imported by the root, so it is in defaultTargets.
Commit: `785c955` (pushed). `lake build` (3496 jobs) exits 0, and `lake test` exits 0.
**Source:** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`, `lem:D`, l.183–186 (the mirror claim), read first-hand. The tie rule is from `def:phi`(3), l.104–105.
It continues the forward half, `slide_eq_switch_product` (`2190619`, note `2026-11-17-lean-slide-switch-C.md`).

## Main theorem (sorry-free)

```lean
theorem reverse_slide_eq_switch_product (D : Finset Cell) (U : List Cell) (h : Cell)
    (hD : IsDiag D) (hS : StdFrom D (U ++ [h])) :
    switchSideRev D U h = jdtSideRev true U h
```

* `switchSideRev D U h = lettersOf (carryRev |U| (chainOf D (U ++ [h])))`. Here `carryRev m` applies `s_m` first, then `s_{m-1}`, …, `s_1`, using the
  unchanged tier-B `switchAt`/`switchL`.
* `jdtSideRev true U h` is defined independently. `jdtStepRev` looks up the west `(r,c-1)` and south `(r-1,c)` neighbours of the hole in the list of `U`-cells
  in **decreasing** entry order, and moves the hole to the larger entry, which is the smaller index. `larger := false` is the lesion. At the boundary, truncated
  subtraction makes the "neighbour" the hole itself. That cell is never in the list, and the docstring says so.
* The hypothesis `StdFrom D (U ++ [h])` says that `U_1 … U_m F_j` is a standard filling over `D`, so `h` is an outer corner of `D ∪ U ∪ {h}`.

## Route: direct mirror (route 2). Chosen in the first minutes
I chose it because C1 `switchL_step` applies unchanged: the triple `D ⊂ D∪{x} ⊂ D∪{x,z}` is the same lemma with the roles renamed. The rotation
route (route 1) would have needed complement-of-diagram bookkeeping and a transport of `jdtStep` under rotation, which is no cheaper. The whole module took about 40 min.

| piece | declaration | content |
|---|---|---|
| switch side | `carryRev_cons` | `carryRev (m+1) (a :: ch) = switchAt (a :: carryRev m ch) 1`, i.e. `s_2…s_m` act on the tail and `s_1` comes last |
| switch side | `carryRev_chainOf` | the leftward carry turns `chainOf D (U ++ [h])` into `chainOf D (carryWRev h U)`, which is again standard. Induction from the left of `U`, so `s_1` compares `U_1` with the final hole |
| jdt side | `jdtRev_cons_irrel` | a cell that is not `≤` the hole is never looked at again (the hole only moves weakly SW) |
| jdt side | `jdtStepRev_adj_head` | an adjacent (W/S) head cell has index 0, i.e. the largest entry, so it is chosen |
| jdt side | `not_ge_rev_of_not_adj` | **where standardness is used**: if the next `U`-cell `y` is not W/S-adjacent then `y ≰ h`. Otherwise the successor of `y` towards `h` lies in the shape and is not `h`, not placed (invariant), and not in `T \ {y}` (a diagram that would then contain `y`) |
| jdt side | `jdtRev_eq_carryV` | for fuel ≥ \|U\|+1, the reverse slide returns `carryV h U.reverse`. Invariant `JdtInvRev`: `StdDown T V`, `insert h (T ∪ P)` is a diagram, `h ∉ T`, and no placed cell is `≤` the hole |
| bridge | `carryV_append`, `carryV_reverse` | the two intermediate words agree |
| setup | `stdFrom_append`, `stdDown_of_stdFrom`, `isDiag_union_of_stdFrom`, `fresh_of_stdFrom` | the initial invariant from `StdFrom D (U ++ [h])` |
| corollary | `reverse_slide_eq_switch_product_le5` | `reverseSlideCheck 5 = true` (162 instances), **derived from the theorem**: the kernel decides only `hypCheckRev 5` |
| corollary | `reverse_slide_eq_switch_product_le6` | `reverseSlideCheck 6 = true` (495 instances), the same way |
| control | `control_rev_both_neighbours` | hole (1,1) with both a W and an S neighbour. The final word is pinned explicitly |
| lesion | `reverse_slide_eq_switch_product_larger_false` | the general statement with `larger := false` is **false**: the control instance meets both hypotheses (`lesion_rev_smaller`) |
| lesion | `reverseSlideCheck_larger_fails` | `reverseSlideCheck 4 false = false`, kernel-decided |

**Two honest corrections made during the session.**
- My first control value was a hand computation, and it was wrong: I forgot that `U_1` slides too, after `U_3` moves. `decide` refuted it, and the
  pinned value was corrected to the one computed by the definitions.
- My first whole-check lesion was at `|ν| ≤ 3`. `decide` showed it was vacuous, because a hole with both a W and an S neighbour first occurs at `|ν| = 4`. It was moved to `|ν| ≤ 4`.

The instance counts 54 / 162 / 495 for `|ν| ≤ 4 / 5 / 6` match the forward enumeration. That is expected: both enumerate pairs (inner diagram, nonempty standard filling).

No mathematics beyond the tex was needed. The paper's "mirror argument" goes through as stated, and there is no gap.

## `#print axioms` (scratch files `/tmp/lc/Ev.lean`, `/tmp/lc/Ax2.lean`, outside the tree, each with a sorry control in the same run)

```
reverse_slide_eq_switch_product              [propext, Classical.choice, Quot.sound]
carryRev_chainOf / jdtRev_eq_carryV          [propext, Classical.choice, Quot.sound]
reverse_slide_eq_switch_product_le5 / _le6   [propext, Classical.choice, Quot.sound]
reverse_slide_eq_switch_product_larger_false [propext, Classical.choice, Quot.sound]
reverseSlideCheck_larger_fails               [propext, Classical.choice, Quot.sound]
control_rev_both_neighbours                  [propext, Classical.choice, Quot.sound]
slide_eq_switch_product (forward, re-checked) [propext, Classical.choice, Quot.sound]
sorry_control                                [sorryAx]
```

The module has 0 sorries.

## Registry
`q585-ybe-train-vs-commutor.json`, node `G4-slide-is-switch-product-formal`: approach and verification now say **forward + reverse**. Its `lean`
field stays the forward name, and the reverse name is recorded in `approach`. `registry_validate.py` and trustcheck (`--files-dir .`) both report
`OK: … is valid`. The parent `H1-G4-rotation-map-is-switching` **stays `proved`**. That needs C3.

## Not done: C3
The dictionary to `GrowthConfluence` (U↔D, F↔R, `initWord p m`, switches = Flips at `L = switchL`, with the forward and reverse products as the
by-rows and by-columns flip orders) has not been started. I did not state it this session.
