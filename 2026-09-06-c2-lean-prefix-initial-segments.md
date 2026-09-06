# LEAN 2026-09-06 c2 — prefix/suffix of an ascending enumeration

**Project:** `projects/lean/tworow_d4_kernel`
**File:** `TworowD4Kernel/PrefixSignSum.lean` (reachable from the root module `TworowD4Kernel`)
**Registry:** `proofs/registry/fock-ribbon-sign-operator.json`
**Paper proof:** `proofs/2026-09-05-Q85-literal-gcd.tex`, Proposition `prop:N`, §"The prefix sign sum"

## Target

Close `Q85-prefix-sign-sum-lean-general-gap` by proving the one lemma the node named:

```lean
theorem take_ascList_toFinset {r : ℕ} (hT : T ⊆ range k) (hr : r ≤ T.card) :
    ((ascList k T).take r).toFinset = S ↔
      S ⊆ T ∧ S.card = r ∧ ∀ x ∈ S, ∀ y ∈ T \ S, x < y
```

## What builds sorry-free

Six new declarations. `lake build` (2978 jobs) and `lake test` both exit 0. Zero sorries in the
file.

| declaration | content |
|---|---|
| `eq_of_lt_sdiff` | two equinumerous subsets of `T`, each strictly below its complement in `T`, are equal |
| `take_ascList_spec` | the length-`r` prefix of `ascList k T` **has** the initial-segment property |
| `take_ascList_toFinset` | **the target**, both directions |
| `eq_of_gt_sdiff` | mirror of `eq_of_lt_sdiff`, order flipped |
| `take_reverse_ascList_spec` | the length-`j` prefix of `(ascList k T).reverse` **has** the final-segment property |
| `take_reverse_ascList_toFinset` | the mirror target |

`#print axioms` on all six:

```
'TworowD4Kernel.eq_of_lt_sdiff' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.take_ascList_spec' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.take_ascList_toFinset' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.eq_of_gt_sdiff' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.take_reverse_ascList_spec' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.take_reverse_ascList_toFinset' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Exactly the standard three. No `native_decide`, no new axiom.

## The shape of the proof

The lemma splits cleanly into *existence* and *uniqueness*, and only the first half touches lists.

**Existence** (`take_ascList_spec`) is four one-liners once the right Mathlib lemmas are named:
`List.take_sublist` gives `Nodup` and containment; `List.toFinset_card_of_nodup` with
`List.length_take` gives the count; and the separation is `List.pairwise_append` applied to
`ascList = take ++ drop` — an element outside the prefix but inside `T` is in the *suffix*, and
`ascList_pairwise_lt` then compares them. No induction on `r` was needed. LEAN.md anticipated one;
the `take ++ drop` decomposition removes it.

**Uniqueness** (`eq_of_lt_sdiff`) has no lists in it at all. If `S ≠ S'` with `#S = #S'`, then
neither contains the other (`Finset.eq_of_subset_of_card_le`), so pick `x ∈ S \ S'` and
`y ∈ S' \ S`; each set's separation hypothesis places one strictly below the other, and `asymm`
closes it. This is the half that makes the characterisation an *iff* rather than a property.

Mathlib has no `take`/`sort` interface of this shape — I checked `List.Sorted`,
`Finset.sort_sorted`, `Finset.orderIsoOfFin`, and grepped for `toFinset_take`. Nothing.

## The correction: the gap named itself wrongly

The node stated the gap as *one* lemma, "after which both cases collapse by
`∑_{B⊆C}(-1)^{|B|}=0`". Reading the paper proof against the Lean definition of `rho` shows that
was wrong. `ρ_T = ascList k T ++ k :: (ascList k D).reverse`, and the paper's case `k ∈ S̄`,
regime `r > m+1`, reads its prefix off the **descending tail** — the `r-m-1` *largest* elements of
`D`. `take_ascList_toFinset` says nothing about that. So I proved the mirror
`take_reverse_ascList_toFinset` as well.

Worth recording as a pattern: the gap statement was written by looking at the *first* case of the
paper proof and assuming the second was symmetric. It is symmetric — but symmetry between two
cases is a reason to prove *two* lemmas, not a reason to prove one.

## The remaining gap, narrower

`prop:N` for general `k` is still open in Lean; `prefixSignSum_eq_three/_four/_five` (`decide`) are
still the only general-`S` statements. But the remaining work is no longer order-theoretic:

1. **Split the prefix across the append.** `((rho k T).take r).toFinset` in the three regimes
   `r ≤ m`, `r = m+1`, `r > m+1`, via `List.take_append_eq_append_take`. The two lemmas above then
   identify each piece. No new mathematics.
2. **Two sum reindexings — this is where the work now is.** Case `k ∉ S̄` reindexes the `T`-sum
   along `T = S̄ ⊔ B`, `B ⊆ (max S̄, k-1]`; case `k ∈ S̄` along `T = (S̄₀ \ Y) ∪ Y'`, `Y' ⊊ Y`.
   Both are `Finset.sum_nbij'`-shaped.
3. **The alternating sum — no work.** `Finset.sum_powerset_neg_one_pow_card_of_nonempty` is already
   in Mathlib (`Mathlib/Data/Nat/Choose/Sum.lean`). LEAN.md was right that this is routine; it is
   in fact free.

Node `Q85-prefix-sign-sum-lean-general-gap` stays `unclassified` — the general theorem is not
formalised, and nothing here is sorried into existence. New sibling node
`Q85-prefix-sign-sum-lean-initial-segments` is `lean-verified` on `take_ascList_toFinset`.

## Tooling facts, confirmed not rediscovered

- `Finset.sort` stayed out of every statement; `ascList` is `(List.range k).filter` throughout.
- No `decide` was added, so the `#guard`-is-not-`decide` trap did not arise.
- Single-file iteration with `LP.sh` ran at **3.5s per check**, against ~10 min for `lake build`.
  The estimate in LEAN.md (~40s) is conservative by an order of magnitude when only one file
  changes and its imports are warm.
