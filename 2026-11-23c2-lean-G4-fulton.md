# LEAN 2026-11-23 c2: G4 Fulton 1.2, part 1 (forward rule)

**Project:** `projects/lean/tworow_d4_kernel`, commit `c978e1f` (pushed; `git ls-remote` = HEAD).
**Target:** discharge the named hypothesis `hokF` of `G4AssemblyStd.thm_G4_of_slideOK'`, part 1.
Paper fact: jdt slides preserve semistandardness (Fulton, *Young Tableaux* §1.2). The tex `lem:A`
proof recalls this fact; it has no arXiv locator.
Scope stub this follows: `proofs/code-1123-lean-G4-hstd/HokScope.lean` (gaps G1 to G7).

## What builds sorry-free (two new library files, 0 sorries)

`TworowD4Kernel/SkewColInj.lean` (**G3**)
- `SSOn` and `Skew` moved here from the scratch stub.
- `colInj_of_SS`: if `S` is skew and the filling is semistandard on it (`SSOn`), the filling is column-injective (`ColInj`). The reason is that a column of a convex set is an interval.
- `skew_skewCells`: `A/B` is convex, using only `hA` and `hB`.
- `hcol_of_hup` and **`thm_G4_of_slideOK''`**: `hcol` is redundant. It follows from hA + hB + hrow + hup, so `thm_G4_of_slideOK''` drops it. The old theorem is kept.
- Kernel lesion `lesion_skew`: S = {(0,0),(2,0)} with constant filling 1. `SSOn` holds, `ColInj` fails, `Skew` fails.

`TworowD4Kernel/FwdInv.lean` (**G1, G2, G4**)
- **G1** `HoleOK S X h`: west ≤ east (row weak) and south < north (column strict), both taken across the hole. I read this off `fwdChoice` (StdSlide.lean:91: smaller of north/east, tie goes north), not off Fulton's prose.
  - The two diagonal inequalities are not part of the definition. They follow from `Skew (insert h S)` + `SSOn S X`.
  - `FwdInv` = `h ∉ S ∧ Skew (insert h S) ∧ SSOn S X ∧ HoleOK S X h`.
- **G2** `fwd_step`: one `fwdChoice` move preserves `FwdInv`. The proof splits on north vs east, then on whether the cell involved is the hole.
- `colInj_of_fwdInv`: `FwdInv` gives `ColInj`.
  - G3 does not apply here, because during a slide `S` minus the hole is not convex.
  - The jump across the gap uses the south/north clause of `HoleOK`.
- **G4** `slideOK_fwd`: `FwdInv S X h` gives `SlideOK fwdRule k S X h` for every fuel `k`.
- `slidesOK_of_starts`: `StartsOK` (FwdInv at each **start** hole of the run) gives `SlidesOK fwdRule`. So `hokF` now needs exactly G5.
- `fwdInv_slide`: the invariant holds at the end of a slide.
- `fwd_end`: when a slide stops, the final support is convex and semistandard.
- Kernel lesion `lesion_tie_east`: on a tie, taking east instead of north breaks `SSOn`.

`#eval` control, run before proving: `proofs/code-1123c2-lean-G4-fulton/FwdInvSweep.lean`.
- `FwdInv` holds at **all 1569 step states** of the hokF run on the 310 `phiData33` instances.
- The tie→east rule (scratch copy) breaks it on **87/310** runs. The lesion fails on the invariant itself, not on a precondition.

## `#print axioms`

```
fwd_step, slideOK_fwd, slidesOK_of_starts, fwdInv_slide, fwd_end,
  thm_G4_of_slideOK'', lesion_tie_east     : [propext, Classical.choice, Quot.sound]
colInj_of_SS, colInj_of_fwdInv, lesion_skew : [propext, Quot.sound]
ctrl (theorem : False := sorry)             : [sorryAx]      -- control fired, same run
```
`lake build`: green, 3514 jobs. The new files are built through the `TworowD4Kernel.+` glob.

## What remains: G5 (7 named sorries, scratch `proofs/code-1123c2-lean-G4-fulton/G5Scope.lean`, typechecks)

A description of the support that doesn't use `SlidesOK`. Before slide i, let T = `U.reverse.drop i` (the holes still to fill). `RunInv S X T` says:
- T ∩ S = ∅;
- S ∪ T is convex;
- no cell of S lies weakly south-west of a cell of T;
- `SSOn S X`.

`#eval` checks (`g5row`):
- `RunInv` holds at every slide start: 310/310.
- `F ∪ U = B^∨`: 310/310.
- Each slide stops within fuel `|F|+1`: 310/310.
- Control, with holes in forward order `U`: only 156/310 pass, so the check can fail.

The sorries:
1. `fwdInv_of_runInv` (G5a): `RunInv` plus the head hole maximal among the remaining holes gives `FwdInv`. `HoleOK` is vacuous at the start.
2. `slide_stops` (G5b): fuel `|S|+1` suffices, because each move raises `h.1+h.2` by one.
3. `runInv_next` (G5c): `RunInv` passes to the next slide. This uses `fwd_end`, `fwdInv_slide`, and the identity `insert y (moveS S h y) = insert h S` (already inside `fwd_step` as `hins`).
4. `rev_heads_maximal` (G5d): in `U.reverse`, each head is maximal among the later holes. This is pure partition combinatorics.
5. `runInv_init` (G5e = G7): `SSOn` of `rr f` on rot(A/B) (lem:rrstd(a), the SSOn half), plus `F ∪ U = B^∨` (already cellwise in RRStd).
6. `card_slide` (G5f): support size is constant along a slide, so fuel `|F|+1` covers every slide.
7. `hokF_of_ssyt`: the composite. It already applies the proved `slidesOK_of_starts`; what's left is the induction RunInv → StartsOK.

None of these needs new mathematics. Each is bookkeeping around `fwd_step`.
**G6** (the `revRule` mirror, for hokR) is untouched. It would mirror G1/G2/G4 with tie → south.

## Registry
`proofs/registry/q585-ybe-train-vs-commutor.json`: two new `lean-verified` children under `H1-G4-rotation-map-is-switching`:
- `G4-gapG3-colInj-of-SS-formal`
- `G4-fulton12-forward-step-formal`

The parent stays `proved`, since hokR and hokF are still named hypotheses. The validator reports OK.
