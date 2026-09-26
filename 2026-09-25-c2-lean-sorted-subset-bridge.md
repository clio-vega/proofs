# LEAN 2026-09-25 c2 — the sorted ↔ subset bridge `J = J'`

**Status: closed, sorry-free.** Target hit; no drift, no new mathematics.

## Target and project

| | |
|---|---|
| Project | `lean/tworow_d4_kernel` (`clio-vega/tworow-d4-kernel`), Lean 4.30.0 / Mathlib `v4.30.0` |
| New file | `TworowD4Kernel/SortedSubsetBridge.lean` (namespace `SortedBridge`) |
| Commits | `cc0573b` (the bridge), `85054bb` (docstring correction in `MConvexExchange`) — both pushed |
| Paper proof | `proofs/2026-09-20-c1-cylindric-M-convexity.tex`, `lem:JJ` / first sentence of the proof of `prop:perm-mconvex` |

Target declaration:

```lean
theorem SortedBridge.insupp_iff_sorted {ℓ : ℕ} {lhat α : ℕ → ℤ} (hl : IsPart ℓ lhat) :
    MConvexExchange.InSupp ℓ lhat α ↔ SortedBridge.InSuppSorted ℓ lhat α
```

`InSupp` (subset form, from `MConvexExchange`) is `α(S) ≤ Λ_{|S|}` for every `S ⊆ [ℓ]`
together with `α([ℓ]) = |λ̂|`; `InSuppSorted` (the paper's form) is `sort(α) ⊴ λ̂`, i.e.
`Dom (sortDesc α ℓ) lhat` together with the same size condition. The **equal-size
convention of dominance** — the content the paper's one-sentence clause elided — is an
explicit conjunct on both sides, so the `λ̂=(2,1)`, `α=(1,1)` counterexample to `⊇`
recorded in the registry cannot arise.

Downstream, and the actual point of the exercise:

```lean
theorem SortedBridge.sorted_exchange {ℓ : ℕ} {lhat α β : ℕ → ℤ} (hl : IsPart ℓ lhat)
    (hα : InSuppSorted ℓ lhat α) (hβ : InSuppSorted ℓ lhat β) {i : ℕ} (hi : β i < α i) :
    ∃ j, α j < β j ∧ InSuppSorted ℓ lhat (ex i j α)
```

— the M-convex exchange axiom for the **paper's own** `J`, not for the subset-form
surrogate. Yesterday's `MConvexExchange.insupp_exchange` only had the surrogate.

## Sorry count

**Zero.** Nothing in the file is sorried, nothing is postponed, no axiom was wanted.

```
'SortedBridge.insupp_iff_sorted'      depends on axioms: [propext, Classical.choice, Quot.sound]
'SortedBridge.sorted_exchange'        depends on axioms: [propext, Classical.choice, Quot.sound]
'SortedBridge.sortDesc_antitone'      depends on axioms: [propext, Classical.choice, Quot.sound]
'SortedBridge.sortDesc_multiset'      depends on axioms: [propext, Classical.choice, Quot.sound]
'SortedBridge.sum_le_sum_take'        depends on axioms: [propext, Classical.choice, Quot.sound]
'SortedBridge.sortDesc_example'       depends on axioms: [propext]
'SortedBridge.sorted_not_unsorted'    depends on axioms: [propext, Classical.choice, Quot.sound]
'SortedBridge.insupp_sorted_witness'  depends on axioms: [propext, Classical.choice, Quot.sound]
```

Exactly the standard three (`sortDesc_example` needs only `propext`: it is `decide`).
`lake build` green across all 2991 targets, and warning-free in the new file.

## The shape of the proof

The brief predicted the real work would be in "sum of the `r` largest", not in the
dominance bookkeeping. That was right. The file factors as:

1. **An explicit sorting permutation of the indices.** `idxSort α ℓ :=
   List.insertionSort (keyRel α) (List.range ℓ)`, where `keyRel α i j := α j ≤ α i` — a
   preorder, *not* antisymmetric (distinct indices may carry equal entries), which is
   why `Multiset.sort` is the wrong tool and `List.insertionSort` (needing only
   `Std.Total` + `IsTrans`) is the right one. `sortDesc α ℓ c` is then the `c`-th entry
   of `α` read along that permutation, and `0` past position `ℓ`.
2. **The crux**, isolated with no `λ̂` anywhere in it:
   ```lean
   lemma sum_le_sum_take (L : List ℤ) (hL : L.Pairwise (· ≥ ·)) (r) (T : Multiset ℤ)
       (hT : T ≤ (L : Multiset ℤ)) (hc : Multiset.card T = r) : T.sum ≤ (L.take r).sum
   ```
   Induction on `L`, splitting on whether the head lies in `T`. The `a ∉ T` branch needs
   `r ≤ L'.length`, and that is not a hypothesis — it falls out of
   `Multiset.card_le_card`, so the awkward case is vacuous rather than requiring
   positivity.
3. **The two halves of `max_{|S|=r} α(S) = Λ`-free sum of the `r` largest**:
   `asum_le_psum_sortDesc` (`≤`, for every `S`, via 2) and `exists_subset_asum_eq`
   (attainment — the witness is `((idxSort α ℓ).take r).toFinset`, which is why the
   permutation was built on *indices* rather than on values; a multiset of values has no
   preferred set of positions to hand back).
4. The bridge itself. `r > ℓ` is handled separately: `psum (sortDesc α ℓ) r` is constant
   past `ℓ`, equals `psum α ℓ = psum λ̂ ℓ`, and `psum λ̂` is monotone because `λ̂ ≥ 0`.

No Mathlib rearrangement lemma was reusable: `Mathlib.Order.Rearrangement` is about
`MonovaryOn` and permutations of a full index set, not about sums over sub-multisets of
a fixed cardinality. Writing the crux directly was cheaper than coercing it.

## Guards (the discipline items)

The brief's standing worry is
`a-definition-transported-into-Lean-is-unfalsifiable-inside-Lean`: if `sort` is *defined*
in Lean, the type checker cannot see a wrong definition. Three guards, in increasing
strength:

* `sortDesc_example` — `sortDesc` *computes* the descending rearrangement of `(1,3,2)`
  as `(3,2,1,0,…)`, by `decide`. A reversed sort direction, or sorting the indices
  instead of the values, fails this.
* `sortDesc_antitone` **and** `sortDesc_multiset` — these are **theorems**, not
  definitional facts, and together they characterise `sortDesc α ℓ` uniquely among
  functions vanishing off `[ℓ]` (antitone + same multiset of values on `[ℓ]`). So a
  wrong `sortDesc` is falsifiable *inside* Lean. This is the guard that was missing in
  the `AffineAdditive` case.
* **Negative control** `sorted_not_unsorted`, proved rather than asserted: with
  `λ̂ = (3,1)`, `ℓ = 2`, `α = (0,4)`, the *unsorted* dominance condition `∀ r, psum α r ≤
  Λ_r` holds together with the size condition, yet `¬ InSupp 2 λ̂ α` (the largest entry
  `4` exceeds `Λ_1 = 3`). So deleting the word `sort` makes the theorem **false** — the
  check can fail.
* Non-vacuity: `insupp_sorted_witness` exhibits an inhabitant of `InSuppSorted` and a
  genuine move produced by `sorted_exchange`.

## Differential check, widened

`proofs/code-q254-lean/mconvex_exchange_check.py` builds both sets from the paper's
definitions, independently of the Lean file. A new `eq` mode runs the `(EQ)` comparison
alone — it is `O(compositions · 2^ℓ)` rather than the `O(|J|²)` of the exchange checks,
so it goes much further:

| | range | result |
|---|---|---|
| before | `\|λ̂\| ≤ 9`, `ℓ ≤ 4` | 0 disagreements / **164** pairs `(λ̂, ℓ)` |
| now | `\|λ̂\| ≤ 16`, `ℓ ≤ 6` | 0 disagreements / **1804** pairs `(λ̂, ℓ)` |

Command: `python3 proofs/code-q254-lean/mconvex_exchange_check.py eq 16 6`; output
`proofs/code-q254-lean/check-eq-d16-ell6.out` (2m25s). `eq 18 7` does not finish inside
a session and was not run — saying so rather than letting the table imply coverage.

## Registry

`proofs/registry/cylindric-M-convexity.json`:

* `sorted-form-equals-subset-form` : `proved` → **`lean-verified`**,
  `lean = SortedBridge.insupp_iff_sorted`.
* `polymatroid-exchange` : `proved` → **`lean-verified`**,
  `lean = SortedBridge.sorted_exchange`. This is the upgrade the session was for: the
  one gap holding the parent back was the subset/sorted identification, and it is closed.
* `polymatroid-exchange-subset-form` : text corrected — it said the identification was
  not formalised, and that is no longer true.

`trustcheck.py … validate … --files-dir .` → `OK … is valid`, exit 0.
`registry_validate.py` prints ~20 `file not found under /home/clio/projects/proofs`
lines and exits 0; these are the **known root-directory defect** in that script (it
resolves `file:` paths under `projects/proofs` rather than `projects`). Confirmed by
`ls` on the four files I touched, not assumed.

## What is still owed — unchanged by this session

**Murota's symmetric axiom (B-EXC)** is *not* formalised. `sorted_exchange` gives the
one-sided exchange: `∃ j, α_j < β_j ∧ α - e_i + e_j ∈ J`. Murota additionally demands
`β + e_i - e_j ∈ J` **for the same `j`**. The two axioms cut out the same class of sets
by a Murota–Shioura theorem that is neither formalised nor used here; the symmetric form
is covered only by brute force (0 failures / 876317 triples). So `lean-verified` on
`polymatroid-exchange` means *the paper's Prop 6.2 as the paper states it* — the paper
quotes the one-sided form as its definition of M-convexity, as do WZZ and Brändén–Huh —
and it does **not** mean Murota's axiom. That caveat is now written at the parent node,
not only at the child; it was reachable only from the child before, which is how a
caveat gets lost on promotion.

Also unchanged: **Rick's re-review of the `a6c83ed` build is still owed.** His endorsement
predates the erratum to the displayed sign in this very proof, and nothing today
discharges it.

## One thing worth keeping

The brief said `J = J'` "is a one-sentence step carrying an implicit *because*", and
`what-is-not-an-assertion-is-not-checked` says to point Lean at exactly those. It cost
one session and came out clean — no defect in the paper's claim this time. That is worth
recording as such: the rule earns its keep by being *cheap when the step is sound*, not
only by the times it finds a false clause. What the formalisation did surface is smaller
and organisational — the honest caveat in `MConvexExchange`'s docstring became a false
sentence the moment the bridge landed, and would have sat there being copied. A caveat
is a claim about the state of the world elsewhere; it rots the same way a provenance
sentence does.
