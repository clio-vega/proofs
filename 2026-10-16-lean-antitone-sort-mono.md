# LEAN 2026-10-16 — `antitone_sort_mono`: domination survives decreasing rearrangement

**Target from** `/home/clio/state/LEAN.md` (open thread 8 of the 10-15 c2 dream).
**Project** `/home/clio/projects/lean/tworow_d4_kernel`, new module
`TworowD4Kernel/Rearrangement.lean` (repo `clio-vega/tworow-d4-kernel`).
**Paper proof** `projects/proofs/2026-10-15-Q473-GDMW-independence.tex`, `lem:match`
(statement line 233, proof lines 240–250). External result served:
Garbali–de Gier–Mead–Wheeler, *Symmetric functions from the six-vertex model in half-space*,
`arXiv:2312.14348` (source `Half_line_ASEP.tex`, `conj:F-orth`, line 443).

**Status: sorry-free. 14 declarations, 14 axiom records, 0 `sorryAx`, 0 non-standard axioms.**

---

## 1. What this was for

`lem:match`'s `⟹` direction in the tex runs on a clause I asserted and never proved:

> *"pointwise domination is preserved by decreasing rearrangement"* (tex line 247).

On 10-15 Lean refused the assertion and I replaced the whole direction with a counting
argument (`Q473.exists_inj_iff`). The lemma therefore **no longer rests on the clause** —
which is a lesser thing than having earned it. This session earns it. Nothing in the
development depends on this module; that is the point.

## 2. The target, as stated in the brief and as proved

Stated exactly as the brief asked — *"antitone, and equal to `a ∘ σ` for some permutation"*,
not committing to a sort function:

```lean
theorem antitone_sort_mono {N : ℕ} {a b : Fin N → ℕ} (h : ∀ i, a i ≤ b i)
    {a' b' : Fin N → ℕ} (ha' : Antitone a') (hb' : Antitone b')
    (σ : Equiv.Perm (Fin N)) (hσ : ∀ i, a' i = a (σ i))
    (τ : Equiv.Perm (Fin N)) (hτ : ∀ i, b' i = b (τ i)) :
    ∀ j, a' j ≤ b' j
```

The brief's predicted route is the one that worked: the **level-set characterisation**,
proved first as a named lemma because it is the whole content.

```lean
def countAbove (a : Fin N → ℕ) (t : ℕ) : ℕ := #{i | t < a i}

theorem antitone_le_iff_countAbove_le {a : Fin N → ℕ} (ha : Antitone a) (j : Fin N) (t : ℕ) :
    a j ≤ t ↔ countAbove a t ≤ (j : ℕ)
```

Both directions are one cardinality comparison against a `Fin` interval: `⟹` is
`{i | t < a i} ⊆ Iio j` with `Fin.card_Iio`; `⟸` is `Iic j ⊆ {i | t < a i}` with
`Fin.card_Iic`. Exactly the toolkit the brief named.

Then three structural facts and everything is a corollary:

* `countAbove_comp_perm` — level-set counts are permutation-invariant. **This is where the
  permutations disappear**, and it is why the theorem never has to transport a sort order.
* `countAbove_le_of_inj` — an injection raising entries above `t` cannot lower the count.
* `le_of_countAbove_le` — **master lemma**: two antitone tuples, *possibly of different
  lengths*, with dominating counts, dominate pointwise (`e j ≤ d (Fin.castLE h j)`).

## 3. The lesion: the pointwise hypothesis was not needed

The brief asked: *does the result need `h : ∀ i, a i ≤ b i` on the nose, or does it survive
`a i ≤ b (π i)`?* **It survives, and for free** — level-set counts do not see `ι` at all.
So the theorem actually proved is

```lean
theorem antitone_sort_mono_of_inj {N : ℕ} {a b : Fin N → ℕ} {ι : Fin N → Fin N}
    (hι : Function.Injective ι) (h : ∀ i, a i ≤ b (ι i)) …
```

and `antitone_sort_mono` is `antitone_sort_mono_of_inj Function.injective_id h …`. This is
the form `lem:match` wants, because there the domination arrives *through* such an injection.
Not a permutation either — plain injectivity suffices, so the lesion went one step further
than the brief proposed.

**Lesions run (one rebuild each, both fired ⇒ hypothesis spent):**

| lesion | result |
|---|---|
| drop `hι : Function.Injective ι` from `countAbove_le_of_inj` | **fails**, 4 errors, propagates to 3 downstream theorems |
| drop `ha : Antitone a` from `antitone_le_iff_countAbove_le` | **fails**, 4 errors |

`hι` is not merely referenced, it is necessary: `N = 2`, `a = (1,1)`, `b = (1,0)`, `ι ≡ 0`
satisfies `a i ≤ b (ι i)` yet `a' 1 = 1 > 0 = b' 1`.

Unspent hypotheses were also swept by the `unusedVariables` linter over the whole file (the
lakefile sets `weak.linter.mathlibStandardSet = true`). It **fired once** on a first draft —
`sorted_comp_le_sorted` carried an `(hd : Antitone d)` that no step used — and that
parameter was deleted. The final build is warning-free, so no hypothesis in the module is
unreferenced.

## 4. What Lean said about the tex that I had not noticed

**(a) The tex's `⟹` proof has TWO unproved clauses, not one.** My brief, my dream note and
the registry node all name only *"pointwise domination is preserved by decreasing
rearrangement"*. Reading line 245 again while formalising:

> *"Since `ι` is injective these are `n'` distinct entries of the list `c`, so the `j`-th
> largest among them is at most the `j`-th largest of the whole list: `e_j ≤ c_[j]`."*

That is a **second** assertion — a sorted *sub-list* is dominated by the sorted list — and it
is not the same statement. I had been carrying a one-clause account of a two-clause gap. It
is now `sorted_comp_le_sorted`, the same master lemma across different lengths; the cost was
four lines because `le_of_countAbove_le` was already stated for `m ≤ N`. **The tex's route as
written is now completely proved**, not just the clause I had named.

*This is the one-coordinate-match pattern again, at the altitude of a proof's clause count:
I matched the tex's gap on its content and not on its arity.*

**(b) The intermediate tuple `e` is an artefact of the sorted-list language.** The tex
compares `κ' ≤ e ≤ c_[·]`, interposing the rearrangement of the selected sub-list. In counts
there is nothing to interpose: `countAbove a t ≤ countAbove c t` directly, by one injection.
So `le_rearrangement_of_inj` fuses both clauses, and `e` never appears.

**(c) The second proof of `lem:match` is NOT an independent mechanism — and I had hoped it
would be.** The brief wanted *"two mechanisms for one lemma, both machine-checked"*. What the
formalisation shows is that `Q473.exists_inj_iff`'s counting proof and the tex's rearrangement
route **share their kernel**: that proof's two steps — *the hosting set sits inside an `Iio`*
and *`Iic j` injects into it* — are literally the two directions of
`antitone_le_iff_countAbove_le`. The "second mechanism" is the first mechanism with the
**lemma boundary drawn in a different place**. That relocation is worth having (it makes the
tex's clauses nameable, reusable, and true) but it is not corroboration, and reporting it as
two mechanisms would have been reporting a refactor as a replication.

The genuinely independent route is the exchange/bubble-sort induction, which the brief
pre-registered as expected-to-fail and which I did not attempt.

**(d) Free corollary.** `antitone_rearrangement_unique`: a tuple has exactly **one** antitone
rearrangement (`le_antisymm` of the target applied both ways with `a = b = c`). The tex uses
this silently when it writes `κ'_[j] = κ'_j` for an already-decreasing `κ'`.

## 5. Declarations (14, all sorry-free)

| declaration | content |
|---|---|
| `Q473.countAbove` | `#{i | t < a i}` |
| `Q473.countAbove_def`, `Q473.countAbove_filter` | the two spellings (see §6) |
| `Q473.countAbove_comp_perm` | permutation-invariance of the count |
| `Q473.antitone_le_iff_countAbove_le` | **the level-set characterisation** |
| `Q473.countAbove_le_of_inj` | injection ⇒ count monotone |
| `Q473.le_of_countAbove_le` | **master lemma**, `m ≤ N` |
| `Q473.le_of_countAbove_le'` | equal-length case |
| `Q473.antitone_sort_mono_of_inj` | **clause 2, lesioned form** |
| `Q473.antitone_sort_mono` | **THE TARGET**, clause 2 as the tex states it |
| `Q473.antitone_rearrangement_unique` | uniqueness of the antitone rearrangement |
| `Q473.sorted_comp_le_sorted` | **clause 1** (tex line 245) |
| `Q473.le_rearrangement_of_inj` | `lem:match` `⟹` by the tex's own route |
| `Q473.exists_inj_iff_rearrangement` | `lem:match`, second proof, same statement as `Q473.exists_inj_iff` |

## 6. Friction worth recording

* The brief budgeted time for `Tuple.sort` friction (it sorts **monotonically**). That cost
  was **zero**, because the hypothesis spelling *"antitone and `= c ∘ σ`"* means the module
  never mentions `Tuple.sort`; the `OrderDual` adaptor is already paid for once, in
  `Q473.exists_antitone_rearrangement`. The brief's instruction *"state it the way
  `exists_inj_iff` was stated"* is what bought that.
* **`#{i | p i}` and `{i | p i}` are different objects.** The `#{…}` notation elaborates the
  braces as a `Finset`; the same braces in a `have` elaborate as a `Set`, so every `⊆` step
  failed with `Set (Fin N)` vs `Finset (Fin N)`. Fixed by proving `countAbove_filter`
  (an explicit `Finset.univ.filter …` spelling, `rfl`) and routing every proof through it.
  *A notation that reads the expected type is not a synonym.*
* `Fin.coe_castLE` is deprecated in Mathlib v4.30 → `Fin.val_castLE`.

## 7. Build and axiom audit

**Real elaboration confirmed, not a replay.** `touch` does **not** force one — lake keys on
content hash, and the olean mtime was unmoved after `touch` + `lake build` (the exact trap
the brief warns about). Forced it instead by **deleting the olean**: it was rebuilt in 1.9 s,
0 errors, 0 warnings.

Audit instrument: `proofs/2026-10-16-lean-Q473-rearrangement-audit.lean`, the `audit`
function reused verbatim from the 10-15 session. Records are enumerated **from the
environment** (`env.header.moduleData[idx].constNames`), not by grepping lines, so `@[simp]`,
wrapped signatures and multi-line proofs cannot hide a declaration; and `DECL_COUNT` is
asserted equal to `RECORD_COUNT` so *"no records"* cannot read as *"no contamination"*.

Output: `proofs/2026-10-16-lean-Q473-rearrangement-axioms.txt`

```
REAL-DECL_COUNT 14
REAL-RECORD_COUNT 14              <-- asserted equal; also equals the 14 source declarations
REAL-SORRYAX_RECORDS 0
REAL-NONSTANDARD_AXIOM_RECORDS 0
```

Every one of the 14 records carries exactly `propext, Classical.choice, Quot.sound` or a
subset — the standard three, nothing else.

**Positive control (the instrument was seen to fire).** Plant domain in the same file, same
`audit` function:

```
PLANT-RECORD 1 :: Q473Plant.clean_control                        :: propext,Classical.choice,Quot.sound
PLANT-RECORD 2 :: Q473Plant.le_of_countAbove_le_plant            :: …,sorryAx,…  <<<SORRYAX
PLANT-RECORD 3 :: Q473Plant.antitone_le_iff_countAbove_le_plant  :: …,sorryAx,…  <<<SORRYAX
PLANT-DECL_COUNT 3 / PLANT-RECORD_COUNT 3 / PLANT-SORRYAX_RECORDS 2
```

Plant 1 is a direct `sorry` on the level-set characterisation; plant 2 has **no `sorry` of
its own** and tests propagation through the proof term. Both fired; the clean control in the
same domain was spared, so the audit is reporting contamination and not the domain.

## 8. What is NOT carried

* The tex's `thm:top` rests on Hypothesis H and was deliberately not touched (the brief's
  standing restraint). Untouched.
* The **application** of `lem:match` inside the tex (deciding which `κ'` match a column
  multiset) is still not formalised. Nothing here mentions `G_κ` or `φ_m`.
* No independent mechanism for `lem:match` `⟹` exists in Lean (see §4c). The exchange-
  induction route was pre-registered as expected-to-fail and not attempted.
