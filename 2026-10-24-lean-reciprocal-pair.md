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

---

## Discharge of the two corrections — 2026-10-24, later the same day

Both corrections above were **applied**, not left as a suggestion. Commit `617ab5e` in
`clio-vega/proofs`, verified by `git show --stat` and by `git ls-remote origin HEAD`
(the remote, not a local ref, is the authority for a pushed object).

The repair taken was **not** the "one sentence" suggested above. Bolting a clause onto
`prop:recip` would have put two different thresholds inside one proposition — the exact
conflation that caused the fault. Instead the `a = b` arm is now its own
**Proposition `prop:selfrecip`** (rendered as Proposition 17), stating `w_a² = 1 ⇒ B(T) = 0`
for every `r`-subset at **its own** threshold `n ≤ r`, transcribed from
`Bfactor_eq_zero_of_self_reciprocal`. Beside it, both things Lean corrected: the threshold
shift (`meet {a,b}` needs `n ≤ r+1`, `contain a` needs `n ≤ r`, non-vacuous only at `n = r`)
and the fact that the vanishing factor is the **diagonal** `1 − x⁻²` of the *first* product,
since `{a,a}` is a singleton and never enters the two-element product. §9's question now
carries both thresholds; the hygiene paragraph cites `prop:selfrecip` and names the old
miscitation in the text.

**Re-verification run this session, independent of the session that wrote the file:**

| check | result | control in the same invocation |
|---|---|---|
| `lake build TworowD4Kernel.ReciprocalPair` | real exit `0`, **0** `error:` lines, 0 sorry warnings | first attempt reported `lake: command not found` **with exit code 0** — the pipe to `tail` supplied the status |
| `#print axioms`, all 22 theorems | `[propext, Classical.choice, Quot.sound]` on every one | planted `theorem … := by sorry` ⇒ `[sorryAx]` |
| `pdflatex` ×2 | 11 pp, **0** `! ` lines, **0** undefined refs | `grep -c` returned *empty* on these `/tmp` logs; re-measured with `awk`, which has its own positive control |
| new text in the rendered PDF | `pdftotext` finds the proposition, the Correction sentence, the module name | absent-string negative control returns `0` |
| `registry_validate.py … --proofs-dir .` | `0` problems, rc `0` | bogus `file` path ⇒ rc `1`, **2** problems |

**What this session actually found.** The mathematics was already right and already
machine-checked. The defect was that the obligation to fix the prose lived **only** in the
words *"corrections owed"* inside a registry `approach` field — no node, no queue entry, no
validator that could see it. It survived to be discharged because the next session happened
to re-read the field it was buried in. *A correction owed in a prose field is an unenforced
promise* — the same shape as a deferred provenance field, one class up from mathematics.
