# LEAN c2 2026-10-23 — Q559's boundary case: the hypothesis I had never written

**Project:** `projects/lean/tworow_d4_kernel`
**Module:** `TworowD4Kernel/BoundaryEmptyIndex.lean` (new) + `TworowD4Kernel/RhoShiftVanishing.lean` (existing, 0 sorry)
**Paper:** `projects/proofs/2026-10-23-Q559-vanishing-converse.tex`, `prop:boundary` (L409–L418),
`cor:criterion` (L392–L406), `thm:alt` (L129–L144), `eq:Psidef` (L112–L118).

## PRE-REGISTRATION — written before the first `lake build` of this module

Timestamp of this block: written at session start, before any Lean in `BoundaryEmptyIndex.lean` compiled.

| # | lesion | PREDICTION | why |
|---|---|---|---|
| L1 | drop `r ≤ M.card` from the converse (`thm:main`) | **IRREMOVABLE** | `prop:boundary` says it is necessary; needs an explicit witness with `M.card < r` and `rhoShift θ` injective |
| L2 | `M.card < r` vs `M.card ≤ r - 1` over `ℕ` | **SPLIT** | equivalent for `r ≥ 1`; at `r = 0`, `n ≤ 0-1` truncates to `n ≤ 0` which is TRUE at `n=0` while `0 < 0` is FALSE |
| L3 | `r = 0` edge | **ABSENT (vacuous)** | `M.card < 0` is unsatisfiable, so the boundary arm never fires at `r = 0`; `Ψ_{(1^0)}` must be the empty product `1`, not `0` |
| L4 | require `M` nonempty | **ABSENT** | should not be needed anywhere; the L1 witness itself is expected to be `M = ∅` |
| {L1,L3} | `r = 0` **with** `r ≤ M.card` dropped | **PAIR PREDICTION: jointly CONSISTENT (i.e. L1 becomes REMOVABLE at `r=0`)** | at `r=0` the card hypothesis is vacuous, so the lesioned statement should still be true — L1's irremovability should be witnessed only at `r ≥ 1` |

Grades are recorded below the build, not here.

## Honest scope (written before claiming anything)

The Lean layer is the **determinant/alternant** layer (`thm:alt`'s index set and summand shape).
`thm:alt` itself — the step from the head-injection definition `eq:Psidef` to the alternant sum —
is **NOT formalised** and is not in scope. `thm:main` (the non-vanishing induction) is **NOT
formalised**; it enters the criterion as an explicit *hypothesis* binder.
So the child node may move; the parent `clio-Q553-vanishing-converse-Q559` **stays `proved`**.
**A Lean-verified summand is not a Lean-verified sum.**

---

## Result: sorry-free, 21 declarations, standard axioms only

`TworowD4Kernel/BoundaryEmptyIndex.lean` (new, imports `RhoShiftVanishing`), built against
`leanprover/lean4:v4.30.0` + Mathlib `v4.30.0`. Full-project `lake build`: 3480 jobs,
`Build completed successfully`; forced re-elaboration of the module alone emits **zero**
warnings, including `unusedVariables`.

### The paper's one sentence is two Lean facts

`prop:boundary`'s proof is one sentence with a parenthesis, and the parenthesis is a
*different theorem*. Both are formalised, because they live over different index sets:

| paper reading | index set | Lean | mechanism |
|---|---|---|---|
| "no injection `h` from `r` slots into an `n`-element set" | `eq:Psidef` | `no_head_injection_of_card_lt`, `sum_headInjections_eq_zero_of_card_lt` | pigeonhole, `Finset.card_le_card_of_injOn` |
| "`M` has no `r`-subset in `thm:alt`" | `thm:alt` | `powersetCard_eq_empty_of_card_lt`, `sum_powersetCard_eq_zero_of_card_lt`, `psiAlt_eq_zero_of_card_lt` | `Finset.powersetCard_eq_empty` |

Both hold for an **arbitrary summand** `f`, so neither needs the summand's actual shape.
That is why the head-injection half is a complete formalisation of `prop:boundary`
*without* `thm:alt`: whatever `Ψ`'s summand is, a sum over an empty index set is `0`.

### Main declarations

* `psiAlt` — the alternant sum of `thm:alt`, with the ordering of each `T` abstracted as a
  parameter `ord` (so the vanishing statements hold for *every* ordering convention).
* `psiAlt_eq_zero_of_card_lt` — **arm 1**, `prop:boundary`. Proof: one `rw`.
* `psiAlt_eq_zero_of_not_injective` — **arm 2**, `thm:vanish` lifted from one determinant
  to the whole sum, over `RhoShiftVanishing.det_zpow_eq_zero_of_not_injective`.
* `psiAlt_eq_zero_of_criterion` — `cor:criterion`'s `⇐` direction as a single proposition,
  with the paper's `j < k` hypothesis verbatim.
* `criterion_iff` — `cor:criterion` in full, with `thm:main` as an explicit hypothesis
  binder `hmain : ∀ θ, r ≤ M.card → Injective (rhoShift θ) → psiAlt … ≠ 0`.
  **The `r ≤ M.card` is in the binder, not in the prose. That is the whole point.**

### `#print axioms`

```
psiAlt_eq_zero_of_card_lt                 : [propext, Classical.choice, Quot.sound]
psiAlt_eq_zero_of_criterion               : [propext, Classical.choice, Quot.sound]
criterion_iff                             : [propext, Classical.choice, Quot.sound]
main_without_card_hypothesis_is_false     : [propext, Classical.choice, Quot.sound]
sum_headInjections_eq_zero_of_card_lt     : [propext, Classical.choice, Quot.sound]
control_criterion_rhs_fails               : [propext, Classical.choice, Quot.sound]
psiAlt_rank_zero                          : [propext, Classical.choice, Quot.sound]
```

No `sorryAx`. Standard three only. **Sorries: 0.**

## Lesion grades — against the predictions above

| # | prediction | GRADE | the Lean declaration that grades it |
|---|---|---|---|
| L1 | IRREMOVABLE | **IRREMOVABLE — confirmed, and by a contradiction, not an absence** | `main_without_card_hypothesis_is_false`: a converse stated *without* `r ≤ M.card` yields `False`. Witness `r = 1`, `M = ∅`, `θ = 0`. |
| L2 | SPLIT | **SPLIT — confirmed** | `card_lt_iff_le_pred` (holds for `r ≥ 1`) + `card_lt_iff_le_pred_fails_at_rank_zero` (fails at `r = 0`: `ℕ`-truncation makes `n ≤ r-1` *true* at `n = r = 0` where `n < r` is false) |
| L3 | ABSENT (vacuous) | **ABSENT — confirmed** | `card_lt_zero_is_unsatisfiable` (the hypothesis cannot hold at `r = 0`, so nothing is claimed) + `psiAlt_rank_zero` (`= B ∅ · (1/V ∅)`, the empty product) + `psiAlt_rank_zero_eq_one` (`= 1 ≠ 0`, so a statement that *did* claim vanishing at `r = 0` would be false) |
| L4 | ABSENT | **ABSENT — confirmed** | `sum_powersetCard_eq_zero_of_card_lt` carries no nonemptiness binder at all; `psiAlt_eq_zero_at_empty` applies it at `M = ∅`, which is where L1's witness lives |
| {L1,L3} | jointly CONSISTENT | **CONSISTENT — confirmed** | `main_without_card_hypothesis_holds_at_rank_zero`: at `r = 0` the lesioned statement is **true**. So L1's irremovability is witnessed only at `r ≥ 1` — at rank `0` the card hypothesis is vacuous and dropping it is harmless. |

Five for five on the first build. That is the pattern I distrust, so:

### Positive controls — because every lesion above would also hold of a constantly-zero `psiAlt`

* `control_psiAlt_rank_one`: at `r = 1`, `M = {5}`, `θ = (3)` the sum evaluates to
  `5^{m_1} = 5^2 = **25**`, computed from the definition, not via any theorem above. So
  `psiAlt` really evaluates its summand and arm 1 is not a theorem about nothing.
* `control_criterion_rhs_fails` (**the separating one**): at that same point *neither*
  disjunct of `cor:criterion` holds (`M.card = 1 = r`, and `rhoShift θ` is injective) and
  the value is nonzero. So `criterion_iff`'s right-hand side is not always true and arm 1's
  hypothesis `M.card < r` carries load rather than being satisfied for free.

### Detector controls — the ones that could have manufactured this green

* **Planted sorry.** Replaced arm 1's proof with `sorry` and rebuilt: `lake build` still
  printed **`Build completed successfully`** (the trap, measured again), while the
  `declaration uses 'sorry'` warning fired and `#print axioms` reported **`sorryAx`**. So
  the two detectors I graded from can in fact fire. Restored; `grep -c sorry` → `0`.
* **Declaration census, both directions.** Source set (21 names, parsed from the file) vs
  Lean's environment set, filtered to the namespace: `comm -23` **empty**; `comm -13` is
  `{headInjections.eq_1, psiAlt.eq_1}` — the two auto-generated equation lemmas of the two
  `def`s, i.e. accounted for, not missing declarations. `21` real declarations.
* **Sorry grep with a live positive control in the same invocation**: `sorry` → 0 while
  `theorem` → 21 in the same read, so the counter was not returning the reassuring value.

## The finding

`psiAlt_eq_zero_of_card_lt` is **one `rw`**. `no_head_injection_of_card_lt` is four lines
of pigeonhole. The hypothesis `n ≥ r` that Q559 was missing for a day was that cheap to
force — and the `424/424` computation that graded the conjecture could not have found it at
any sample size, because its generator drew `n ∈ r…r+2`. The evidence was a restatement of
the unwritten hypothesis.

What Lean adds beyond "the proof is short": **L1's grade is a contradiction, not a gap.**
Stated without `r ≤ M.card`, the converse is not unproved — it is *refutable*, and the
refutation is four tokens of witness. A binder Lean can refute is a binder prose can lose.

And the pair {L1, L3} is the honest qualifier the paper does not carry: `n ≥ r` is
necessary **for `r ≥ 1`**. At `r = 0` it is vacuous, and the lesioned statement is true. So
the sharp form of `prop:boundary`'s necessity claim has a rank hypothesis of its own.

## Registry

`proofs/registry/rick-q488-converse.json`, node
`clio-Q559-boundary-n-ge-r-was-omitted`: `proved` → **`lean-verified`**, with the four
declaration names in `lean` and `lean_file` set.
`python3 code/registry_validate.py proofs/registry/rick-q488-converse.json --proofs-dir .`
→ `OK: … is valid`, **0 problem lines** (and it printed a positive `OK:` line, so it ran).

The parent `clio-Q553-vanishing-converse-Q559` **stays `proved`**: `thm:alt` and `thm:main`
are not formalised. A Lean-verified summand is not a Lean-verified sum.
