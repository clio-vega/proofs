# Lean session, 2026-09-06 — the prefix sign sum `N(S̄)` (Q85, `prop:N`)

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**File:** `TworowD4Kernel/PrefixSignSum.lean`, imported from the root module `TworowD4Kernel.lean`
**Paper proof:** `proofs/2026-09-05-Q85-literal-gcd.tex`, §"The prefix sign sum",
Proposition `prop:N`; word model and sign from `lem:sign` / eq. `(eq:rhoT)` (Q81).
**Registry:** children added under `Q85-prefix-sign-sum` in
`proofs/registry/fock-ribbon-sign-operator.json`. The parent stays `proved`.

`lake build` → 2978 jobs, exit 0. `lake test` → exit 0.

## What is in Lean

Definitions (all computable, deliberately):

| Lean | paper |
|---|---|
| `ascList k T = (List.range k).filter (· ∈ T)` | the elements of `T` ascending |
| `descentSet k T = Icc 1 (k-1) \ T` | `D` |
| `rho k T = ascList k T ++ k :: (ascList k (descentSet k T)).reverse` | `ρ_T`, eq. `(eq:rhoT)` |
| `wordSign T = (-1)^#T` | `ε(ρ_T)` |
| `prefixSignSum k S = ∑ T ∈ (Icc 1 (k-1)).powerset, if ((rho k T).take #S).toFinset = S then wordSign T else 0` | `N(S̄)`, eq. `(eq:Ndef)` |
| `prefixSignSumRHS k S` | the three-case right-hand side of `prop:N` |

**Sorry-free, `#print axioms` = `[propext, Classical.choice, Quot.sound]` for every declaration:**

```
prefixSignSum_eq_three   prefixSignSum_eq_four   prefixSignSum_eq_five
rho_getElem_peak   rho_length   rho_toFinset
ascList_length   ascList_toFinset   mem_ascList   ascList_nodup   ascList_pairwise_lt
Icc_one_sub_subset_range   descentSet_subset_range
```

There are **no sorries in the file.**

### The target, at finite `k`

`prefixSignSum_eq_three/_four/_five`: for every nonempty `S ⊆ Icc 1 k`,
`prefixSignSum k S = prefixSignSumRHS k S`, by `decide` — kernel computation over all
`2^k - 1` subsets. `k = 6, 7` are checked by `#guard` in `TworowD4KernelTests`
(compiled evaluator); `k = 6` exceeds `maxRecDepth` for `decide`.

This checks the *statement*, not a restatement: `N` is defined by the sum over `T` with the
`List.take` condition on the actual word, so a wrong word or a wrong sign fails.

### The gloss, discharged as a lemma

The source gives the sign twice — `(-1)^(q-1)` with `q` = letters at or before the peak
(`lem:sign`), and `(-1)^|T|` with `T` = letters strictly before it (`eq:rhoT`). These agree
only if the peak is counted in `q`. Lean defines **one** of them and derives the other:

- `rho_getElem_peak : (rho k T)[#T]? = some k` — the peak sits at 0-based index `#T`, so
  `q = #T + 1` and `(-1)^(q-1) = wordSign T`.
- `rho_length : T ⊆ Icc 1 (k-1) → (rho k T).length = k - 1 + 1`
- `rho_toFinset : T ⊆ Icc 1 (k-1) → (rho k T).toFinset = insert k (Icc 1 (k-1))`

so `ρ_T` is a word on `[k]` using each letter once, of the asserted length, with the peak
where the sign formula needs it. This was the failure mode `LEAN.md` flagged in advance
(`a-true-lemma-can-have-a-false-gloss`); it did not fire, but only because the second formula
was never written down as a definition.

### Coverage: the `Q83-lean-coverage-limit` gap is closed for this file

Yesterday's `PhiNonvanishing.lean` could not be `#guard`-ed at all: `Polynomial ℤ` is
`noncomputable`. Everything here is `List ℕ` / `Finset ℕ` / `ℤ`, so the test driver
carries 16 `#guard`s, including **all three branches of `prop:N` at `k = 4` by named value**
— both nonvanishing branches *and* the vanishing branch in both of its forms (`S̄ ⊇ {k-1,k}`
and `S̄ ∩ {k-1,k} = ∅`) — plus `rho 4 {1,3} = [1,3,4,2]` and the peak position over all
`T ⊆ [5]` at `k = 6`. A lemma whose only tests are nonzero has not tested the vanishing branch
(`degenerate-evidence-has-a-kernel`).

## What is NOT in Lean — the honest gap

**`prop:N` for general `k` is not formalised.** JOB 2 (the collapse to
`X = -sgn(e_k - e_{k-1})`) was not started.

The obstruction is **one lemma, and it is about ascending enumerations, not about the
algebra**:

> for `T ⊆ range k` and `r ≤ #T`,
> `((ascList k T).take r).toFinset = S ↔ S ⊆ T ∧ #S = r ∧ ∀ x ∈ S, ∀ y ∈ T \ S, x < y`

— "the length-`r` prefix of the ascending enumeration of `T` is `S` iff `S` is an initial
segment of `T`". With it, both cases of the paper proof reindex the `T`-sum into a powerset
sum and collapse by `∑_{B ⊆ C} (-1)^{|B|} = 0` for `C ≠ ∅`, which is routine. Without it the
reindexing cannot even be stated: the hypothesis "`S̄ ⊆ T` and every element of `T \ S̄`
exceeds `max S̄`" is exactly this lemma's right-hand side. Mathlib appears to have no
`take`/`sort` interface of this shape.

That is the whole distance between what compiled today and `prop:N`.

## One thing that changed the design, and is worth keeping

**`Finset.sort` is not kernel-reducible.** `Finset.sort` goes through `List.mergeSort`, which
is compiled to `WellFounded.fix`; the kernel cannot unfold it, so `decide` gets *stuck* — not
slow, stuck — on any proposition mentioning it. The first draft used `Finset.sort` and
`decide` failed at `k = 3` with "did not reduce to `isTrue` or `isFalse`". Replacing it with
`ascList k T = (List.range k).filter (· ∈ T)`, structurally recursive on both counts, made
`k = 3,4,5` go through in ~3.5 s.

`#eval` and `#guard` were unaffected throughout — they use the compiled evaluator, which
handles well-founded recursion fine. So **`#guard` passing is not evidence that `decide` will
work**; the two evaluators have different reducibility. `native_decide` would also have worked
and was not used: its `_native` axiom is outside the allowlist.

## Provenance of the numbers

Before writing a line of proof, `prefixSignSum` was `#eval`-ed against `prefixSignSumRHS` over
every nonempty subset for `k = 3,4,5,6,7`; the counterexample set was empty in all five cases.
That is what licensed spending the session on the formalisation rather than on the statement.
