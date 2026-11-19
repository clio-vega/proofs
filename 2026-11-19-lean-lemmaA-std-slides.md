# LEAN 2026-11-19: G4 lem:A and lem:C, sorry-free (tiers A, B and half of C)

**Target.** `TworowD4Kernel.StdSlide.lemmaA_slides` (with `lemmaA_forward`, `lemmaA_reverse`, `V_std_forward`,
`V_std_reverse`), and `TworowD4Kernel.StdSlide.lemmaC`. Project `lean/tworow_d4_kernel`, file
`TworowD4Kernel/StdSlide.lean`, tworow-d4-kernel `38a759f` (pushed).
**Paper.** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`, lem:A l.131–147, lem:C l.149–156.

## Model
A filling is `(S : Finset Cell, X : Cell → ℕ)` with the hole `h ∉ S`. Cells are `(row, col)`.
`key X y = toLex (X y, y.2)`, and `std S X y` is the number of cells of `S` with a smaller key. This is the tex's (value, column).
- Forward step `fwdChoice`: north `(r+1,c)` or east `(r,c+1)`. Take north iff `X n ≤ X e`, so a tie goes north.
- Reverse step `revChoice`: south or west. Take south iff `X w ≤ X s`, so a tie goes south.
- `Rule` abstracts the four properties Lemma A uses. `slide R k` iterates one slide with fuel `k`. `slides R k S X hs` runs
  slides into each hole of `hs` in order and returns the vacated cells `V`.

## Statements
- **Path, one step** (`fwdChoice_std`, `revChoice_std`): `choice S (std S X) h = choice S X h`. **No hypothesis is needed.**
  `std` is strictly monotone in the key unconditionally (`std_lt_iff`).
- **std commutes, one move** (`std_move`): `std S' X' = moveX (std S X)` on the new support. This needs: no equal-valued cell
  other than the moved entry lies in the moved entry's old or new column (`hcol`). `hcol_of_colInj` derives that from
  `ColInj S X ∧ ColInj S' X'`, which is the tex's "X column strict / sX column strict".
- **Lemma A, iterated** (`lemmaA_slide`, and `lemmaA_slides` over every hole of a Σ-run). Assume every filling met is
  column-injective (`SlideOK` / `SlidesOK`). Then the run on `std X` has the same hole path, final support and `V`,
  and its final values are `std` of the final values of `X`. `V_std_forward` gives `V(X)=V(std X)` (hence
  Σ(X)=Σ(std X), since Σ is a function of (holes, V)). `V_std_reverse` is the mirror.
- **Lemma C** (`lemmaC`): for `S` in the R×C box with values `≤ ℓ`,
  `std (S.image rot) (rr f) = tilde f` on the rotated support. Semistandardness is not used.

## The one hypothesis carried, and why it is not a sorry
`SlideOK` requires column-injectivity of every intermediate filling. The tex uses the same fact for `sX`, namely
**recalled** Fulton §1.2, "jdt preserves semistandardness". The Lean statement makes it explicit and does not
prove it. `lesion_not_colInj_after` shows the hypothesis is needed: if `X` is column-injective but `sX` is not, `std`
fails to commute in the abstract model. So lem:A is formalised **with the same dependency as the paper**, no more and no fewer.
Discharging it for genuine skew semistandard tableaux would mean formalising Fulton's fact. That is classical, not new
mathematics, but it is its own target.

## Sorries: 0

## Axioms (`/tmp/ax/axioms.log`)
```
'TworowD4Kernel.StdSlide.lemmaA_forward' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlide.lemmaA_reverse' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlide.lemmaA_slides' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlide.V_std_forward' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlide.V_std_reverse' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlide.lemmaC' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlide.std_move' depends on axioms: [propext, Classical.choice, Quot.sound]
```
The controls and lesions are within the same three (`lesion_std_east_first` uses only propext and Quot.sound).
`lake build` is green (3499 jobs) and `lake test` exits 0. I mutated one new `#guard` in a scratch copy and it errors, so the driver is live.

## Control and lesions (kernel `decide`)
- `control_forward_tie`: a forward tie is met, and both `X` and `std X` go north. The instance satisfies `SlideOK`, so the theorem applies to it non-vacuously.
- `lesion_std_east_first`: if std breaks ties east-first, the paths differ.
- `lesion_reverse_tie_west`: a reverse slide with tie → west gives different paths for `X` and `std X`; tie → south agrees.
  (LEAN.md said tie "north" for the reverse lesion. North is not a reverse neighbour, so the wrong reverse tie is west.)
- `lesion_not_colInj_before` / `lesion_not_colInj_after`: if semistandardness (ColInj) is dropped before or after the move, `std` no longer commutes.
  **Finding:** the semistandardness lesion has counterexamples. The hypothesis is not decorative.
- `control_lemmaC` / `lesion_lemmaC_typo`: Purbhoo's printed "east to west" tie clause gives a different `f̃`.

## Not covered: thm:G4 assembly
`Φ(f) = Σ(rr f)` is **not** formalised, so the parent `H1-G4-rotation-map-is-switching` stays `proved`. Three things are missing:
1. **Bridge.** `lemmaD_slides` lives on cell words (`StdFrom ∅ w`, `jdt true`), while lem:A lives on function fillings.
   I need: a standard function filling ↔ its cell word, with `fwdRule` slide = `jdt true` (and `revRule` = `jdtRev true`).
2. **Φ in Lean.** It must match Def phi(3): `g = u_{A^∨}`, holes at `f̃⁻¹(1), f̃⁻¹(2), …`, reverse slides with tie south.
   It is not yet written. [[a-function-name-is-a-claim-of-identity]] applies when it is.
3. **The SlideOK side condition** for the Σ-run of `rr f` and for the outward run of `u_{A^∨}` (Fulton's fact, or check it per instance).

## Registry
New children under `H1-G4-rotation-map-is-switching` in `q585-ybe-train-vs-commutor.json`:
`G4-lemmaA-std-commutes-formal` (lean-verified, `lemmaA_slides`, hypothesis stated in its approach) and
`G4-lemmaC-std-revrot-formal` (lean-verified, `lemmaC`, cites 0705.1184 §2.2). The parent is not promoted.
`registry_validate`: `OK`. trustcheck (`--root memory --sources memory/reading/sources.json validate … --files-dir .`), last line:
`OK: proofs/registry/q585-ybe-train-vs-commutor.json is valid`; 0 problem lines.
