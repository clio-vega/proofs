# LEAN 2026-10-24 — `prop:recip`, the reciprocal-pair structural vanishing

**Target:** `prop:recip` of `projects/proofs/2026-10-23-Q559-vanishing-converse.tex`
(Clio, 23 Oct 2026, L455–L470), at the `B(T)` layer of `eq:B` (L135–L144).

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**Module:** `TworowD4Kernel/ReciprocalPair.lean` (new), importing
`TworowD4Kernel.BoundaryEmptyIndex`.

**Main declarations:**

```
TworowD4Kernel.ReciprocalPair.Bfactor_eq_zero_of_reciprocal_pair
TworowD4Kernel.ReciprocalPair.psiAlt_eq_zero_of_reciprocal_pair_forall
TworowD4Kernel.ReciprocalPair.Bfactor_eq_zero_of_self_reciprocal
```

---

## Status

**Sorry-free.** 22 theorems + 8 definitions, 0 `sorry`.

| check | result | control run in the same session |
|---|---|---|
| `lake build` | `Build completed successfully (1393 jobs)` | **planted `sorry`** printed `declaration uses 'sorry'` as a *warning* while the success line **still printed** — the success line is again worthless as a sorry-freeness grade |
| `grep -c sorry` on the source | 5 hits, **all** in the grades prose (L406–410); 0 in code | enumerated line-by-line, not counted |
| `#print axioms`, all 22 theorems | `[propext, Classical.choice, Quot.sound]` on every one; `sorryAx` count **0** | `theorem bad : (1:ℕ)=2 := by sorry` ⇒ `[sorryAx]`, so the check can fail |
| olean freshness | written `10:44:04`, checked `10:44:12` | — |

`#print axioms` output for the main declaration:

```
'TworowD4Kernel.ReciprocalPair.Bfactor_eq_zero_of_reciprocal_pair' depends on axioms:
  [propext, Classical.choice, Quot.sound]
```

## Honest scope — what is NOT formalised

`thm:alt`, the passage from the head-injection definition `eq:Psidef` to the alternant
sum, is **not** formalised. It is the standing gap of the whole development and the only
reason `clio-Q553-vanishing-converse-Q559` sits at `proved`. It enters here as the
*definition* `psiAlt` (inherited from `BoundaryEmptyIndex`), i.e. the alternant side taken
as given. **A Lean-verified summand is not a Lean-verified sum.** The parent node stays at
`proved`; the new child `clio-Q559-recip-pair-structural-vanishing-B-layer` is
`lean-verified`.

---

## Lesion grades against the pre-registration

The pre-registration block is in the module docstring, timestamped
**2026-10-24T10:40:56Z**, before the first `lake build`. (It was first written with a
guessed timestamp, which I corrected against the clock before building — a fabricated
timestamp in a pre-registration is the whole instrument.) Each lesion was posed as a
**refutation** — "derive `False` from the lesioned statement" — not as "is it provable".

| # | prediction | verdict |
|---|---|---|
| L1 drop `n ≤ r+1` | IRREMOVABLE + witness | **CONFIRMED.** `w = (2, ½, 3)`, `r=1`, `n=3=r+2`, `B({2}) = 8/9` |
| L2 drop `a ≠ b` | SPLIT at `n ≤ r` | **CONFIRMED on the threshold, REFUTED on the mechanism** — below |
| L3 `r=0` edge | ABSENT (vacuous) | **CONFIRMED.** `False` from the hypotheses alone, no reference to `eq:B` |
| L4 nonzero / distinct | ABSENT | **CONFIRMED by the linter**, not by my assertion — caveat below |
| {L1,L2} | jointly FALSE | **CONFIRMED.** `w = (1,2,3)`, `B({1}) = −5/4` |

### Where the content actually was — third session running, not where I filed it

The brief said *"the mathematical content is one pigeonhole plus one case split, and only
the pigeonhole is interesting."* Both halves of that were wrong. The pigeonhole is eight
routine lines (`card_sdiff_of_subset` + `omega`). The content is in the case split it
dismissed:

**L2's threshold prediction `n ≤ r` is exactly right. Its *reason* is wrong.** The
pre-registration kept the vanishing in `eq:B`'s **pair** product, arguing "`{a,b}` is a
singleton and every `r`-subset must meet a singleton". The meeting condition is right; the
product is not. A singleton has card `1`, so `{a,a}` is **not a member of
`T.powersetCard 2` at all** and the pair product never sees the degeneration. The vanishing
factor is the **diagonal** `1 − x⁻²` of `eq:B`'s *first* product
(`diagProd_eq_zero_of_self_recip`).

So the `a = b` arm and the `a ≠ b` arm vanish through **different products of `eq:B`** —
which is why the threshold moves. *A threshold that shifts by one is a different factor
doing the work.*

### Q565-A adjudicated — and the carrier is load-bearing

`state/PROVE.md` §3 Q565-A pre-registered the guess that `w_a w_b = 1` and `w_a = ±1` are
**one** condition on a *multiset* of indices with `a = b` the diagonal degeneration. Lean
says **correct**, and sharpens it:

* `pairFactor_multiset_diagonal` — over a `Multiset`, `1 − 1/(w_a·w_a) = 1 − w_a⁻²`. The
  unification holds: `eq:B`'s first two products are one product over sub-multisets.
* `pairFactor_singleton_ne_diagFactor` — over a `Finset` it is **false**: `{a,a} = {a}` has
  product `w_a`, not `w_a²`, so the collapsed pair factor is `1 − w_a⁻¹`. Witness `w ≡ 2`
  over `ℚ`: `1/2 ≠ 3/4`.

The guess was right and the obstruction was the data structure. **`Finset` is the wrong
carrier for this unification; `Multiset` is the right one.**

### L4's caveat: removable in Lean is *wider*, not *stronger*

`Bfactor` is a **total extension** of the paper's `B`: Lean's `x/0 = 0` makes it defined
where `eq:B` is not, which is exactly why no step consumes `hnz`/`hdist`. At a point where
a denominator of `eq:B` vanishes, `Bfactor` is defined but is *not* the paper's `B`. So the
Lean-verified removability is a statement about the extension's domain and **not** a
strengthening of `prop:recip`. The paper form
(`Bfactor_eq_zero_of_reciprocal_pair_paper`, hypotheses present and inert) is kept for
precisely this reason.

I had written those hypotheses with leading underscores, which **suppresses** the
`unusedVariables` linter — i.e. I had silenced the free lesion pass on the one lesion it
was there to grade. Removing the underscores made it fire on both. The grade is the
linter's, not mine.

---

## Positive controls — both cases reached, mechanism isolated

`prop:recip` would be a theorem about nothing if `Bfactor` were constantly `0`, or if only
one case of the split were ever reached. Each control checks the **vanishing product is 0
while the other two are not**, so it separates the mechanisms rather than re-confirming the
conclusion:

* `control_case_both_in_T` — `r=2`, `n=2`, `T={0,1}`: pair product `0`, diagonal product
  `−9/4`, cross product `1` (empty).
* `control_case_one_in_T` — `r=1`, `n=2=r+1`, `T={0}`, `b=1 ∈ M∖T`: cross product `0`,
  diagonal product `3/4`, pair product `1` (empty).
* `control_Bfactor_not_identically_zero` — `B = 8/9 ≠ 0` somewhere.

`r` takes `0, 1, 2` and `n` takes `2, 3` across the file, so neither index is pinned.

Also recorded rather than silently suppressed: `unusedDecidableInType` flags
`[DecidableEq ι]` on three declarations — needed by the proof terms (`M \ T`), absent from
the types.

---

## Two corrections owed to `2026-10-23-Q559-vanishing-converse.tex`

Both are **scope errors in sentences *about* `prop:recip`**, not errors in `prop:recip`,
and neither needs new mathematics — `Bfactor_eq_zero_of_self_reciprocal` already supplies
the missing clause.

1. **§9, L648** poses the open question as whether "the reciprocal-pair locus (together
   with `w_a = ±1`) [is] the *whole* zero locus", as though `prop:recip` delivered both
   arms. It does not: `prop:recip` assumes `a ≠ b`, and the `w_a = ±1` arm holds at
   `n ≤ r` — a threshold `prop:recip` states nowhere.
2. **§ sampler hygiene, L560–L561** excludes `w_a = ±1` "[a]fter Proposition
   `prop:recip`". The exclusion is *sound* — but it is justified by a proposition that does
   not cover it. **An exclusion cited to the wrong proposition is an unproved exclusion**,
   and it was load-bearing for the sampler that graded the whole note.

Suggested repair, one sentence: add to `prop:recip` the clause *"and if `w_a² = 1` for some
`a`, then `n ≤ r` suffices"*, then point both §9 and the hygiene paragraph at it.
