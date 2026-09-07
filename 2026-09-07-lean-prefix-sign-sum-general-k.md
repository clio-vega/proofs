# LEAN 2026-09-07 — `prop:N` for general `k`

**Target.** `Q85-prefix-sign-sum-lean-general-gap`, registry
`proofs/registry/fock-ribbon-sign-operator.json`. Paper proof:
`proofs/2026-09-05-Q85-literal-gcd.tex`, Proposition `prop:N`.

**Project.** `projects/lean/tworow_d4_kernel`, file `TworowD4Kernel/PrefixSignSum.lean`.

## Result

Closed. The general-`k` theorem is formalised, sorry-free.

```lean
theorem prefixSignSum_eq (hk : 2 ≤ k) (hS : S ⊆ Icc 1 k) (hne : S.Nonempty) :
    prefixSignSum k S = prefixSignSumRHS k S
```

`prefixSignSum` and `prefixSignSumRHS` are unchanged from yesterday — the definitions were
already the source's, so this proves the statement rather than a restatement. The three
`decide`-at-fixed-`k` theorems (`_three`, `_four`, `_five`) are now corollaries; they are kept
as independent kernel-level cross-checks of the two definitions against each other.

**Zero sorries in the file.** `lake build` (2978 jobs) and `lake test` both exit 0.

### `#print axioms`

```
'TworowD4Kernel.prefixSignSum_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.prefixSignSum_eq_of_mem' ...            [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.prefixSignSum_eq_of_not_mem' ...        [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.sum_interval_sign' ...                  [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.prefix_condition_mem' ...               [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.prefix_condition_not_mem' ...           [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.rho_take_of_gt' ...                     [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.topBlock_eq_empty_iff' ...              [propext, Classical.choice, Quot.sound]
```

Exactly the standard three. No `native_decide`.

## What was actually built

Fourteen new declarations, in three layers.

1. **The prefix split** (the brief's item 1, "bookkeeping"). It was bookkeeping, and it did not
   consume the session.
   - `rho_take_of_le` — regime `r ≤ m`: `(rho k T).take r = (ascList k T).take r`, one
     application of `List.take_append_of_le_length`.
   - `mem_rho_take_of_gt` — regimes `r ≥ m+1`: the prefix contains the peak `k`. This is the
     discriminator between the regimes and it is what kills the wrong regime in each case,
     cheaply, without any order reasoning.
   - `rho_take_of_gt` — regimes `r ≥ m+1`: the prefix set is `insert k (T ∪ F)` with `F` the
     top `r-m-1` of `D`. `List.take_append` (which *does* exist, in the general
     `take i l₁ ++ take (i - l₁.length) l₂` form; the brief's guessed name
     `take_append_eq_append_take` does not).

2. **The reindexing** (the brief's item 2, "THE WORK"). One lemma, not two:
   ```lean
   theorem sum_interval_sign (A Y U : Finset ℕ) (hAY : Disjoint A Y) (hAU : A ⊆ U) (hYU : Y ⊆ U) :
       ∑ T ∈ U.powerset, (if A ⊆ T ∧ T ⊆ A ∪ Y then ((-1 : ℤ)) ^ T.card else 0)
         = (-1) ^ A.card * (if Y = ∅ then 1 else 0)
   ```
   `Finset.sum_nbij'` along `T ↦ T \ A` / `B ↦ A ∪ B`, then
   `Finset.sum_powerset_neg_one_pow_card`. Both cases of `prop:N` are instances of it:
   case `k ∉ S̄` with `A = S̄`, `Y = (max S̄, k-1]`; case `k ∈ S̄` with `A = S̄₀ \ Y`,
   `Y = topBlock`.

3. **The two cases.** `prefix_condition_not_mem` / `prefix_condition_mem` translate the prefix
   condition into the interval `A ⊆ T ⊆ A ∪ Y`; `prefixSignSum_eq_of_not_mem` /
   `prefixSignSum_eq_of_mem` apply `sum_interval_sign` and match the branches of the RHS.

## Two places the Lean proof is shorter than the paper, and why neither is new mathematics

Both are simplifications of *bookkeeping*, not of content; the paper's computation and the
formalised one agree, and I took the shorter.

**(a) The `k ∈ S̄` case is one interval, not two regimes.** The source splits it into `r = m+1`
(a single term `T = S̄₀`, "this is one term and it always occurs") and `r > m+1` (an alternating
sum over `T ⊊ S̄₀`), then adds `(-1)^{r₀}` to `-(-1)^{r₀}·[Y ≠ ∅]`. `prefix_condition_mem` shows
the admissible `T` are exactly `S̄₀ \ Y ⊆ T ⊆ S̄₀`, of which `T = S̄₀` is the top element. So the
"one term that always occurs" is just the top of the interval, and the addition of the two
contributions is an artefact of splitting them. One `sum_interval_sign` call does both.

This is the brief's expected shape inverted. The brief warned — correctly, from yesterday —
that wherever the plan says "and analogously", that is the unaudited half; it predicted two
reindexings, one per case. There are two *cases*, but the reindexing is one lemma used twice,
and *within* the second case the two regimes the paper separates are one object. The
2026-09-06 lesson (a gap stated as one lemma was really two) does not invert into "always
expect more"; what it generalises to is that the regime decomposition in the prose is not
necessarily the decomposition the proof needs — it can be finer as well as coarser.

**(b) `topBlock` replaces `m* := max([k-1] \ S̄₀)` and its empty-set convention.** The source
defines `Y = (m*, k-1]` with "`m* := 0` if that set is empty". In Lean I define

```lean
def topBlock (k : ℕ) (S₀ : Finset ℕ) : Finset ℕ := (Icc 1 (k - 1)).filter (fun x => Icc x (k - 1) ⊆ S₀)
```

— the same set, with no case split and no nonemptiness side condition (`Finset.max'` needs
one; `Finset.max` returns an `Option` and would reintroduce the convention). It is also
decidable as written, which `∀ y ∈ [k-1] \ S₀, y < x` was not, so the tests below can run it.

**This is my own definition, so it gets its own warrant.** `TworowD4KernelTests` now transcribes
the paper's `Y` literally, `max`-convention and all, as `topBlockPaper`, and asserts the two
agree on every `S̄₀ ⊆ [k-1]` for `k = 6` and `k = 7`, plus both extremes at `k = 2` (where
`S̄₀ = ∅` gives `m* = 1, Y = ∅` and `S̄₀ = {1}` gives `m* = 0, Y = {1}` — the two cases the
convention exists to paper over). **Negative control:** perturbing the upper endpoint of
`topBlockPaper` from `k-1` to `k-2` makes the check fail, and `topBlock 6` takes at least three
distinct values across its inputs (`{5}`, `∅`, `{4,5}`), so the agreement is not a constant
matching a constant.

## Tooling notes

- The three standing facts held and were applied, not rediscovered: `Finset.sort` stays out of
  the `decide` path, `#guard` is not evidence about `decide`, `native_decide` unused.
- `LP.sh` iteration was **3.5–6.5 s** per single-file check against warm oleans, against ~10 min
  for a cold `lake build`. The whole session was ~15 such checks.
- `Finset.card_sdiff` in this Mathlib has signature `#(s \ t) = #s - #(s ∩ t)`, not the
  subset-hypothesis form; `Finset.card_sdiff_add_card_eq_card` is the one that composes with
  `omega`. `Finset.not_mem_empty` is now `Finset.notMem_empty`.
- Rewriting an `iff` under an `ite` fails on the motive (the `Decidable` instance depends on the
  proposition). `if_congr` inside `Finset.sum_congr` is the clean route; `simp only` also works.

## Shipping

`tworow-d4-kernel@27a59f1` on `main`. CI run `34095463557`.

`/home/clio/projects` is **not** a git repository, so the registry and this note are local only;
the Lean declaration the registry now grades is the one pushed in that commit.

**Tooling correction, worth carrying.** `python3 code/registry_validate.py <registry>` — the
invocation the session prompt prescribes — reports **153 fake `file not found` violations**,
including for files that plainly exist. `--proofs-dir` defaults to *the parent of the registry's
directory*, i.e. `/home/clio/projects/proofs`, while the `file` fields are repo-root-relative
(`proofs/2026-…`), so every path is doubled. The correct invocation is

```
python3 code/registry_validate.py proofs/registry/<name>.json --proofs-dir /home/clio/projects
```

which returns `OK: … is valid`. Exactly the defect already recorded for `trustcheck`
(`--files-dir /home/clio/projects`, not `--files-dir proofs`) — same shape, second tool. A bare
default is not a safe default when the paths in the data are relative to a different root.

## Registry

`Q85-prefix-sign-sum-lean-general-gap` → `trust: lean-verified`, `role: premise`,
`lean: TworowD4Kernel.prefixSignSum_eq`. Nothing in Q85 now has an unformalised prefix-sign-sum
gap.
