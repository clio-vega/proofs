# LEAN 2026-11-22: G4 assembly, links L1, L3 and L4 sorry-free (Tier A and Tier B); L2 and L6 remain

**Target.** Brief `state/LEAN.md` (11-22). This session discharges the links of the `thm_G4_of_slideOK` assembly stub.
Project `lean/tworow_d4_kernel` (github clio-vega/tworow-d4-kernel). It adds two modules, `TworowD4Kernel/ChainInj.lean` and `TworowD4Kernel/RankList.lean`, and imports both from the root.
Lean commits are `4b3e201` (L1) and `99539b0` (L3+L4). Both are pushed, and `ls-remote origin main` = `99539b0`.
**Paper.** `proofs/2026-11-13c2-Q585-G4-evacuation.tex` (`lem:C`, `lem:D`, `def:phi`(1)); Purbhoo arXiv:0705.1184 §5.1 Cor 5.3.

## Tier A: L1, no sorries (`ChainInj.lean`)
- `chainOf_inj` is a pure list lemma. Let `w₁` have no repeated cell and avoid `κ`, and let `w₂` be any word with `phi (chainOf κ w₁) = phi (chainOf κ w₂)`. Then `w₁ = w₂`.
  Only `w₁` needs a hypothesis. The induction is on `w₁`. Index 1 of the chain gives `insert x κ = insert y κ`, and `x ∉ κ` then forces `x = y`. Lengths are not assumed: if one word is longer, index `|short|+1` compares `∅` with a nonempty set.
- `stdFrom_rowsW` says `rowsW` preserves `StdFrom`. It is the `hstd'` step inside `rows_aux'`'s induction, now stated as its own lemma. It is what makes the output of `rowsW` repetition-free.
- **L1** `rowsW_eq_colsW`: for `StdFrom ∅ w` with `|w| = m+p`, `rowsW p m w = colsW m 0 p w`. The proof is `lemmaD_slides`, then `stdFrom_rowsW` and `stdFrom_nodup`, then `chainOf_inj`.
- **Lesion (kernel):** `[a,b,a]` and `[a,b,b]` have the same chain (`lesion_chainOf_not_inj`). So `chainOf_inj_needs_nodup` holds: dropping `Nodup` makes the statement false.

## Tier B: the tilde-rank lemma, then L3 and L4, no sorries (`RankList.lean`)
**Statement decided before writing it (as the brief asked).** `holesWith t` is `(range n).filterMap (j ↦ rc.find? (t · = j))` with `rc = rot(cs)`.
It enumerates `rot(A/B)` exactly when `t` is injective on `rot(A/B)`. Its range always lies in `[0, |A/B|)`, and injectivity together with that range gives surjectivity by counting.
For `t = tilde`, injectivity is the same as **column-strictness of `f` on `A/B`**. In Lean that is `ColInj`, carried to `rr f` on the rotated support and then to `std` (injective by `std_lt_iff`), with `lemmaC` turning `tilde` into `std (rr f)`.
- **Consumed from `lemmaC`:** only its output `std (S.image rot) (rr ℓ f) = tilde` on `S.image rot`. Its hypotheses are InBox and values `≤ ℓ`; the second holds by `le_foldr_max`, since `ℓ` is Φ's `foldr max 0`.
- **Two new named hypotheses**, both standing hypotheses of `def:phi`(1) and not new mathematics: `hbox : InBox R C (A/B)` and `hcol : ColInj (A/B) f` (`f` is semistandard, so it is strictly increasing down columns).
  **Finding:** the 11-21 stub had neither, and L3 is **false** without `hcol`. See the lesion below. `hstd` does not imply it: `holesWith` is always nodup, so `StdFrom` says nothing about how many holes there are.
- `rankList_spec` (generic): let `r` be injective on `rc` with values `< |rc.toFinset| = n`. Then `rankList rc r n` has `toFinset = rc.toFinset` and `idxOf y = r y` for every `y ∈ rc`.
- `holesWith_spec`, then **L3** `L3_holes_toFinset` (`F.toFinset = rot(A/B)`) and **L4** `L4_std_eq_valOf` (`std F.toFinset (rr f) = valOf F` on `F`).
  The helper lemmas are `skewCells_nodup'` (for general `B`), `colInj_rr`, `std_injOn` and `std_lt_card`.
- **Lesion (kernel):** `lesion_not_colInj` takes `A=(1,1)`, `B=∅`, the 2×2 box and `f ≡ 1`. InBox holds, but `holesWith` has length 1 rather than 2, and the L3 equation fails.

## Assembly stub (scratch, `proofs/code-1122-lean-G4-assembly/AssemblyStub.lean`)
The stub now has `hbox` and `hcol`. **L0, L1, L3, L4 and L5 are discharged** with library lemmas. It type-checks, and its one `declaration uses sorry` warning comes from 2 remaining sorries:
- **L2 (glue)**: Φ = `(U'.toFinset, gU U ∘ valOf U')` with `U' = (colsW …).drop p`. It should follow from `phi_values_destd`, `slides_colsW_take` and `slides_fuel_eq` (fuel R+C ↔ |U|+1, on the whole triple). **Not attempted this session** because the hour ran out.
- **L6 + final combination**: from L1 and the two window forms, derive `V.reverse = U'`, then the labels `row(U[m−1−k])+1 = gU U (valOf U' V[k])`. Not attempted.
The named `SlidesOK` hypotheses (`hokR`, `hokF`) stay named. The stub is **not** in the library and is **not** registered.

## Build, axioms
`lake build` is green (3508 jobs) and `lake test` exits 0. `#print axioms` output is in `proofs/code-1122-lean-G4-assembly/axioms.log`. Of its 10 lines, 9 equal `[propext, Classical.choice, Quot.sound]`; `chainOf_inj` gives `[propext, Quot.sound]`, a subset.
Main declarations:
```
'TworowD4Kernel.ChainInj.rowsW_eq_colsW' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.RankList.holesWith_spec' depends on axioms: [propext, Classical.choice, Quot.sound]
```

## Registry
There are two new lean-verified children under `H1-G4-rotation-map-is-switching`: `G4-L1-chain-injective-rowsW-eq-colsW` (`lean` = `TworowD4Kernel.ChainInj.rowsW_eq_colsW`)
and `G4-L3L4-tilde-rank-holes-enumerate` (`lean` = `TworowD4Kernel.RankList.holesWith_spec`). The parent stays **proved**, since L2, L6 and the SlidesOK discharges remain.
`registry_validate`: `OK`. trustcheck (from `memory/`, `--deployment ../code/clio.json --sources reading/sources.json`, `--files-dir /home/clio/projects`) ends with `OK: … is valid (status: in-progress, deployment: clio)`.

## Next LEAN
L2 (glue only) and then L6. If L2 needs more than the three named lemmas, record that as a finding.
