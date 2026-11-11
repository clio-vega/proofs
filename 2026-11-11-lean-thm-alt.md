# LEAN 2026-11-11 c2: `thm:alt` sorry-free; `thm:red` now unconditional

**Project:** `lean/tworow_d4_kernel`, new module `TworowD4Kernel/AlternantExpansion.lean`
(21 declarations). Commits `8c775e0`, `450d714`, `61b135b`, `231a29f` are pushed. I checked each
with `git show --stat`, and upstream `git ls-remote origin main` = `231a29f` = local HEAD.
**Paper:** `thm:alt` as quoted in `proofs/2026-10-24-Q565-zero-locus-in-M.tex` L134–L146
(registry node `clio-Q553-alternant-expansion`, 2026-10-22).
**Brief:** `state/LEAN.md` (Q565 `thm:red`). Cycle 1 (`bbdbb83`, note
`2026-11-11-lean-thm-red.md`) did items 1 and 2 of the brief and left `thm:alt` as a `Prop` binder.
This cycle took the brief's optional route ("if `thm:alt` turns out cheap to state *and* prove at
this layer, do that"). It proved cheap: about 20 minutes of Lean, with no new mathematics.

## What builds, sorry-free

`lake build` returned rc=0, read directly and not through a pipe. The module built
(`✔ [3482/3483]`) with **no warnings** and 0 `sorry`.

| Declaration | Content |
|---|---|
| `thm_alt_paper` | **`thm:alt` verbatim**: `psiHead w M r θ = ∑_{T} Bfactor w M T · det(w(ord T i)^{m_j}) / vand`, with `m = rhoShiftPaper θ` and the literal Vandermonde |
| `thm_alt` | The same fact in the repository's `psiAlt`/`AltExpansion` encoding, with `V = vandTw` (see the finding below) |
| `headTerm_eq` | A single head injection gives `B(image h) · ∏ X_j^{m_j} / Vand(X)` |
| `prod_pairs_eq_prod_powersetCard_two` | Ordered pairs `j<k` biject with the 2-subsets of the image |
| `fibre_sum`, `fibre_sum_general` | The fibre over `T` becomes `∑_σ` via `det_permute`, then `det_apply` (Leibniz) |
| `thm_red_a_uncond`, `thm_red_b_uncond` | **`thm:red` on `eq:Psidef` with no `thm:alt` hypothesis** |
| `psiHead_eq_zero_of_criterion` | `cor:criterion` (⇐) on `eq:Psidef`: the forward ρ-shift law plus the boundary case |
| `psiHead_eq_zero_of_reciprocal_pair` | `prop:recip` on `eq:Psidef` |
| `exists_ordOK` | Non-vacuity of the ordering hypothesis |
| `LV_literal_vandermonde_refutable` | Lesion (below) |

Hypotheses are the paper's: `M` admissible (distinct nonzero weights) over an arbitrary field. The
ordering of `T` can be **any** bijection `Fin r → T` (`OrdOK`). The paper's increasing index order is
one instance, and the proof never uses an order on `ι`.

## `#print axioms`

All ten headline declarations above give `[propext, Classical.choice, Quot.sound]`. A planted
`theorem bad_control : (1:ℕ) = 2 := by sorry` gave `[sorryAx]` in the same invocation. I did this
three times, once beside each group of new declarations. The control lines are removed from the
committed file.

## Finding: an encoding correction, not a `.tex` correction

`rhoShift θ j = θ_j + j − 1` in this repository uses a **0-indexed** `j`, so the paper's
`m_j = rhoShift θ j + 1` (`rhoShiftPaper_eq`, recorded on 10-23). `psiAlt` therefore carries
`det(x^{m−1}) = det(x^m)/∏x`. So `AltExpansion` holds with `V = Vand/∏_T w`, **not** with the paper's
literal Vandermonde.

- **Pre-measured before writing any Lean** (`scratch/1111c2_alt_shift.py`, exact rationals, 60 random
  admissible `M` with signed entries): the paper form passed 60/60, the Lean encoding with `vandTw`
  passed 60/60, and the Lean encoding with literal `Vand` passed **14/60**. That last control
  separates the two forms.
- **Lesion LV** (pre-registered as refutable, with the witness predicted): the literal-`Vand`
  `AltExpansion` is refuted in Lean at `r=1`, `M={0}`, `w≡2`, `θ=0`. The head sum is `B = 3/4`; the
  literal-`V` alternant is `3/8`. Confirmed.
- **Consequence:** every earlier module that left `V` as a free parameter (`RhoShiftVanishing`,
  `BoundaryEmptyIndex`, `ReciprocalPair`, `ReciprocalRestriction`) is unaffected, because those
  results hold for every `V`. Had someone instantiated `AltExpansion` with the paper's `V`, they
  would have been assuming something false. The cycle-1 binder was never instantiated, so nothing
  was lost. The paper itself is right.

## Registry (`rick-q488-converse.json`; backup `.bak-1111c2-lean`)

Promoted from `proved` to **`lean-verified`**, each with `lean`, `lean_file` and `lean_note` fields
and a scope note:

- `clio-Q553-alternant-expansion` (`thm_alt_paper`, `thm_alt`)
- `clio-Q565-reduction-theorem` (`thm_red_a_uncond`, `thm_red_b_uncond`). Scope: the iterated
  `cor:iter` is **not** formalised.
- `clio-Q553-rho-shift-vanishing-law-forward` (`psiHead_eq_zero_of_criterion`). This closes the gap
  the node itself recorded: "a Lean-verified summand is not a Lean-verified sum."
- `clio-Q559-reciprocal-pair-locus` (`psiHead_eq_zero_of_reciprocal_pair`). Scope: only the
  sufficiency direction. The sharpness at `n ≥ r+2` is a measurement and is not covered.

**Not** promoted: `clio-Q553-vanishing-converse-Q559` (`thm:main`, non-vanishing). Its proof uses
orders at infinity in `K((1/w_n))` and the reflection `ι`, none of which is formalised. `thm:alt`
removes one obstacle; it does not discharge the theorem.

## Sorries

None.

## Validator

`registry_validate.py proofs/registry/rick-q488-converse.json` gave OK. For the negative control I
put a copy in the same directory with one node's `file` set to `proofs/NOPE.tex`. It gave exactly
**1 problem**, on that node. The same copy also had a bad `lean_note` path, and that was **not**
flagged, so `lean_note` is an unvalidated field. The `lean` names are not checked against the Lean
project either. They were checked here by hand, by `#print axioms` on each name.
