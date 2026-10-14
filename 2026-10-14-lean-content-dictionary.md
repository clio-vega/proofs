# LEAN 2026-10-14 c2 — `prop:dict`, content half: `boxCount_eq_one` is now a theorem

**Target.** Turn `TworowD4Kernel.TableauReduction.Tableau.boxCount_eq_one` from a hypothesis
field into a derived fact.

**Project.** `lean/tworow_d4_kernel`, new module `TworowD4Kernel/ContentDictionary.lean`
(imported by the root `TworowD4Kernel.lean`).

**Paper proof.** `prop:dict` of `proofs/2026-10-10-c2-Q410-korff-cylindric-weight.tex`
(statement at line 210; `(eq:mvec)`, `(eq:korffrel)` inside it; `(eq:Psi)` at line 253).
External: Korff **arXiv:1110.6356** Def. 5.5 / Lemma 5.4; Huh-Kim-Krattenthaler-Okada
**arXiv:2301.13117** Prop. 2.4 and Prop. 2.7. Both IDs carried in the module docstring.

---

## §1 — The split, made before any Lean was written

The brief required `prop:dict` be split at its implication arrows and each half named
*before* starting, because the 10-13 session reported `1 -> 0` on a bundled obstruction. The
split is **three ways**, not two:

| | statement | status |
|---|---|---|
| **D1** | content => box count: `Sum_i theta^(a)_i = alpha_a` for a row-strict filling | **in Lean, sorry-free** |
| **D2** | the path: `x^(b)` = row lengths of entries `< b`, and `x^(b+1) = x^(b) + theta^(b)` | **in Lean, sorry-free** |
| **D3** | the region: `x^(b)` in `R = {x_1 >= ... >= x_h >= x_1 - w}` | **not in Lean — hypothesis** |

The brief predicted a **two**-way split (transpose dictionary / Korff's Lemma 5.4 relation)
and guessed the Lemma 5.4 relation would be the tractable half. Both halves of that guess
were wrong, in an informative way:

* The relation `(eq:korffrel)` `m_i(x+theta) - m_i(x) = theta_i - theta_{i+1}` was **not
  needed at all** for the field. It is a statement about `mult`, and `boxCount_eq_one` never
  mentions `mult`. What the field needed was the *content* sentence, which the paper states
  in passing inside its HKKO Prop. 2.7 paraphrase ("the `a`-th step being a `0/1` vector with
  `alpha_a` ones") rather than in `prop:dict`'s own display.
* The residue is **D3**, the region — geometry, not index bookkeeping.

So the honest report is **`1 -> 1`, with the assumption changed, not removed**: `Tableau.mem`
was already a field, and `toTableau` takes exactly that same hypothesis (`hmem`) and nothing
else. No new assumption is introduced; `boxCount_eq_one` stops being one.

## §2 — Domain check, printed side by side (the brief's hard requirement)

| | paper's `prop:dict` | Lean `Filling` |
|---|---|---|
| rows | `h` rows of a cylindric `lambda` in `Par(h,w)` | `ZMod h`, arbitrary `len : ZMod h -> N` |
| rows strictly increasing | yes | **yes** (`Filling.rowStrict`) |
| columns weakly increasing | yes | **no** |
| cylindric `T(i,j) = T(i+k,j-n)` | yes | **no** |
| `lambda` in `Par(h,w)` | yes | **no** |
| entries bounded by `n` | yes | **no** |

**The domains differ, and that difference is the session's finding.** `Filling` is strictly
larger, so D1 and D2 are *stronger* than the paper's claim. Row-strictness is the **only**
tableau condition the content half consumes — visible as the single place it is spent, in
`rowCount_le_one`. The three dropped conditions are consumed entirely by D3.

This is not a convenience. The paper presents the content statement as part of the cylindric
package; it is not — it is a statement about any row-strict filling of any shape. The
cylindricity was carried along because the *English* sentence travelled inside a paragraph
about cylindric tableaux. Same mechanism as the "strip" and "the index map" findings of the
last two sessions: a noun need not declare its domain, a typed binder must.

## §3 — What builds, and the numbers

**Declarations: 26.** All sorry-free.

D1: `Filling` (structure), `rowCount`, `rowCount_le_one`, `boxes`, `content`,
`content_eq_sum_rowCount`, `occ`, `boxCount_occ`, `boxCount_occ_eq_one`.
D2: `rowLen`, `xOf`, `rowLen_succ`, `xOf_succ`.
Constructor + composite: `toTableau`, `rootMultiplicity_Psi_of_filling`.
Non-vacuity: `witF`, `witF_content`, `witF_mem`, `witT`, `witT_occ`, `witT_nuT`,
`witT_natJ`, `witT_multNat`, `witT_Psi`, `witT_rootMultiplicity`.
Ablation: `ablation_rowStrict`.

The headline statements:

```lean
theorem Filling.boxCount_occ (a : N) : boxCount (F.occ a) = F.content a
theorem Filling.xOf_succ (b : N) : F.xOf (b + 1) = shift (F.xOf b) (toNatWord (F.occ b))

def toTableau (F : Filling h) (w : Z) (hc : forall a : Fin n, F.content a.1 = 1)
    (hmem : forall a : Fin n, Mem h w (F.xOf (a.1 + 1))) : Tableau h n w

theorem rootMultiplicity_Psi_of_filling {R} [CommRing R] [IsDomain R] [CharZero R]
    (hh : 2 <= h) (F : Filling h) (w : Z) (hc : forall a : Fin n, F.content a.1 = 1)
    (hmem : forall a : Fin n, Mem h w (F.xOf (a.1 + 1))) :
    rootMultiplicity (1 : R) ((toTableau F w hc hmem).Psi R) = n
```

**Sorries: 0.** Whole-project build: `lake build` from
`/home/clio/projects/lean/tworow_d4_kernel`, **Build completed successfully (3464 jobs)**,
`LAKE_EXIT=0`, **0 error lines**, **0** ``declaration uses `sorry` `` lines, and 0 warnings
naming `ContentDictionary.lean`.

**`#print axioms`: 27 records, 0 non-standard.** 26 declarations of this module plus
`TableauReduction.Tableau.rootMultiplicity_Psi_eq_nuT` (the parent chain it feeds). Every
record within `{propext, Classical.choice, Quot.sound}`; `Filling` itself needs `propext`
alone. Evidence: `proofs/2026-10-14-lean-content-dictionary-axioms.txt`.

> **Record count 27, line count 35.** Four records wrap across three lines each. A line count
> reads 35 and is wrong. Both numbers are printed because that is the §4 control; the records
> are parsed by the `'<name>' depends on axioms:` pattern, not by line.

## §4 — Controls

**Plant A — fabricated declaration names. 2 armed, 2 fired.**
`...ContentDictionary.Filling.boxCount_occ_eq_content` (deliberately inside this file's own
family — a plausible rename of the real `boxCount_occ`) and
`...ContentDictionary.rootMultiplicity_Psi_of_tableau` (a plausible sibling of the composite).
Both: `error(lean.unknownIdentifier): Unknown constant`. `EXIT=1`.

**Plant B — a `sorry` at the dependency root (`rowCount_le_one`, where `rowStrict` is spent).
Contaminated 12 of 18 checked; spared 6.**

Contaminated: `rowCount_le_one`, `boxCount_occ`, `boxCount_occ_eq_one`, `xOf_succ`,
`toTableau`, `rootMultiplicity_Psi_of_filling`, and **all six `witT_*` declarations**.
Spared: `content_eq_sum_rowCount`, `rowLen_succ`, `witF_content`, `witF_mem`,
`ablation_rowStrict`, `TableauReduction.Tableau.rootMultiplicity_Psi_eq_nuT`.

**The sparing is the informative half, and here it reports against me.** On 10-13 the plant
spared `wit_nuT` and `wit_Psi` — kernel computations corroborating by a second mechanism.
This time it reached *every* `witT_*` declaration, including `witT_Psi` and `witT_multNat`.
Reason: `witT := toTableau witF ...`, so the non-vacuity controls route through the chain and
are **not** independent corroboration of it. What survives as a second mechanism is
`witF_content` and `witF_mem` — the kernel deciding the content and the region membership of
the raw filling, with no reference to D1 — together with the pre-existing
`TableauReduction.wit_Psi`, which this plant also spares. A `Tableau`-level non-vacuity
control cannot be independent of a theorem that builds the `Tableau`; to stay independent it
has to live on the raw data, which is what `witF_content` does.

**Ablation — row-strictness is load-bearing, and dropping it makes D1 false, not unproved.**
`ablation_rowStrict` exhibits `h = 2`, row `0` of length `2` with both boxes carrying entry
`0`: the row-sum of row counts is `2` while the box count of the Boolean word is `1`. The two
expressions in the statement are literally the two sides `boxCount_occ` equates, with
`rowCount` unfolded, so it ablates the theorem and not a paraphrase.

**A reading-method fault, caught in this session.** My first check for the plant's sorry was
`grep -c "declaration uses 'sorry'"` -> **0**, on a file that demonstrably contains one. Lean
4.30 writes the warning with **backticks**: ``declaration uses `sorry` ``. A wrong pattern
returns `0`, and `0` is exactly what a clean file returns — the failure is indistinguishable
from success. The build-log count in §3 uses the corrected pattern
`declaration uses .sorry.`.

## §5 — Scope: what Lean carries vs what the nodes assert

* `korff-hkko-dictionary` (trust `proved`) asserts `prop:dict` entire: the two parameter
  identifications `n_K = h`, `k_K = w`, both routes (shape parameters via HKKO Prop. 2.4, and
  the state count), `(eq:mvec)`'s non-negativity and `Sum_i m_i = w`, and `(eq:korffrel)`.
  **Lean carries none of those five.** It carries D1 and D2, which are the part of the
  dictionary the *Lean chain* consumes, and which `prop:dict`'s own display does not state.
  The parent stays at `proved`.
* `cyclic-ascent-lemma` asserts 5 things; the brief measured Lean as carrying 3, and this
  session did not change what that node asserts. It stays at `proved`.
* `q410-no-cyclic-psi-weight` asserts 4 things, Lean carries 1. **Declined for promotion a
  seventh consecutive time**, on the same count that made the first six right.

Per §3 of the brief, the grade goes to a **scope-stated child**:
`content-dictionary-formalised`, `lean-verified`, under `korff-hkko-dictionary`.

## §6 — What did not close, with its exact statement

**D3.** For `lambda` in `Par(h,w)` and `T` in `CRST(lambda;h,w)` a cylindric row-strict
tableau, the vector `x^(b)_i = #{j : T(i,j) < b}` lies in
`R = {x in Z^h : x_1 >= x_2 >= ... >= x_h >= x_1 - w}` for every `b`. In Lean:
`forall a : Fin n, Mem h w (F.xOf (a.1 + 1))`, the `hmem` argument of `toTableau`. This is
HKKO **arXiv:2301.13117, Prop. 2.7** proper, and it is where the column condition, the
cylindric identification and `lambda in Par(h,w)` are spent. Formalising it needs a model
that carries the column condition and the wrap, which `Filling` deliberately does not.

It was already an assumption of `Tableau` before this session, so the chain is no weaker;
what changed is that the single remaining assumption is now **geometric**, with a named
external source, rather than a bundled combinatorial one.
