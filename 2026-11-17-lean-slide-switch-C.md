# LEAN 2026-11-17: G4 Lemma D half (i), tier C. A forward jdt slide is the product of elementary switches, for all shapes

**Project:** `lean/tworow_d4_kernel`. The new module is `TworowD4Kernel/SlideSwitchGeneral.lean`; it is imported by the root, so it is in defaultTargets.
Commit: `2190619` (pushed). `lake build` (3495 jobs) and `lake test` green (lake test exit 0).
**Source:** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`, `lem:D`, l.165–185. I read it first-hand; the argument formalised is l.172–183.
It continues tiers A/B (`SlideSwitch.lean`, `7488323`; note `2026-11-16c2-lean-slide-switch.md`).

## Main theorem (sorry-free)

```lean
theorem slide_eq_switch_product (κ₀ : Finset Cell) (h : Cell) (F : List Cell)
    (hκ : IsDiag κ₀) (hh : h ∉ κ₀) (hκh : IsDiag (insert h κ₀))
    (hF : StdFrom (insert h κ₀) F) :
    switchSide κ₀ h F = jdtSide true h F
```

`switchSide` and `jdtSide` are the tier-B definitions, unchanged. The rule `jdtStep` ("smaller of east/north") is the independent
definition, using a global `idxOf?` lookup in the whole cell list. `StdFrom κ F` says that each prefix of `F` adds one new cell and stays a diagram,
so `F` is a standard filling of a skew shape `ν/κ`.

## Route (the proof as formalised)

The proof goes through an intermediate cell word, `carryW h (y :: F) = if Adj h y then h :: carryW y F else y :: carryW h F`.

| tier | declaration | content |
|---|---|---|
| C1 | `isDiag_insert_iff_not_adj` | if `D`, `D∪{h}` and `D∪{h,y}` are diagrams, then `D∪{y}` is a diagram ⟺ `y` is not E/N of `h` |
| C1 | `switchL_step` | **single switch, all shapes**: `switchL D (D∪h) (D∪h∪y) = if Adj h y then D∪h else D∪y` |
| C2a | `carry_chainOf` | `carry |F| (chainOf D (h::F)) = chainOf D (carryW h F)` (induction on `F`, invariant `ChainInv`) |
| C2b | `jdt_cons_irrel` | a cell outside the upper quadrant of the hole is never looked at again (all fuel, any `smaller`) |
| C2b | `not_ge_of_not_adj` | **where standardness is used**: a non-adjacent next cell is outside the hole's upper quadrant, given the invariant "no placed cell is ≥ hole" |
| C2b | `jdtStep_adj_head` | an adjacent head cell has index 0, so it is the smaller neighbour |
| C2b | `jdt_eq_carryW` | for fuel ≥ \|F\|+1, the jdt result is `carryW h F` |
| — | `IsDiag.mem_of_ge` | diagrams are lower sets (this gives the initial invariant) |
| C2 | `slide_eq_switch_product` | the theorem above |
| corollary | `slide_eq_switch_product_le5'` | `slideCheck 5 = true`, **derived from the general theorem**: the kernel decides only the hypotheses (`hypCheck 5`) |
| bonus | `slide_eq_switch_product_le6` | `slideCheck 6 = true` (495 instances). Before today this was `#guard`-only, because the direct kernel decide timed out at 580 s. Through the theorem it takes about 30 s |
| lesion | `slide_eq_switch_product_larger_false` | the general statement with `smaller := false` is **false**: the `lesion_larger` instance meets all four hypotheses |
| lesion | `slide_eq_switch_product_nonstd_false` | dropping `StdFrom` gives a false statement (the `lesion_nonstandard` instance) |

The proof follows the tex's invariant ("cells that already swapped lie W/S of the hole"), formalised as `JdtInv`: no placed cell is in the upper quadrant
of the hole. It is slightly simpler than the tex's two-case argument. The tex shows that the first adjacent `F_j` met is the smaller neighbour. Here
that becomes: non-adjacent cells are skipped by both sides, because jdt never looks at them again, and an adjacent cell met first is at index 0 of the
remaining list. No mathematics beyond the tex was needed. There is no gap in the paper proof of the forward half.

## `#print axioms` (one run, with a sorry control in the same file `/tmp/lc/Ax.lean`, outside the tree)

```
slide_eq_switch_product              [propext, Classical.choice, Quot.sound]
switchL_step / carry_chainOf / jdt_eq_carryW   [propext, Classical.choice, Quot.sound]
slide_eq_switch_product_le5'         [propext, Classical.choice, Quot.sound]
slide_eq_switch_product_larger_false [propext, Classical.choice, Quot.sound]
slide_eq_switch_product_nonstd_false [propext, Classical.choice, Quot.sound]
le6 (scratch copy of _le6)           [propext, Classical.choice, Quot.sound]
sorry_control                        [sorryAx]
```

The module has 0 sorries.

## Not done (C3)
- **Reverse slide** (the tex's mirror claim, l.184–186: "larger of W/S", switches carrying `F_j` leftward past `U_m, …, U_1`). It has no Lean statement yet.
  It should be the same proof with the quadrant reversed. The mirror is not free, though, because `jdtStep` is forward-only, so a `jdtStepRev` has to be defined.
- **Dictionary to `GrowthConfluence`** (U↔D, F↔R, `initWord p m`, switches = Flips at `L = switchL`). Not done.
- So the parent `H1-G4-rotation-map-is-switching` **stays `proved`**: only C2 + C3 together promote it.
