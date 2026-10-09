# For Robin — 2026-10-09 LEAN: two-row two-part SSYT, `∃!`, sorry-free

Short version: `lem:mono`'s *combinatorial* half is now machine-checked, both halves
(uniqueness **and** existence), 13 declarations, 0 sorry, standard three axioms.

**Lean:** `TworowD4Kernel/TwoRowZeroOne.lean` —
https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/TwoRowZeroOne.lean
**Report:** https://github.com/clio-vega/proofs/blob/main/2026-10-09-lean-two-row-zero-one-unique.md

```lean
theorem SemistandardYoungTableau.existsUnique_le_one_card_ones_eq
    (μ : YoungDiagram) (k : ℕ) (h2 : μ.rowLen 2 = 0)
    (hk1 : μ.rowLen 1 ≤ k) (hk2 : k ≤ μ.rowLen 0) :
    ∃! T : SemistandardYoungTableau μ,
      (∀ i j, (i, j) ∈ μ → T i j ≤ 1) ∧ #{c ∈ μ.cells | T c.1 c.2 = 1} = k
```

Three things worth your time:

1. **The scope boundary is real and I have not crossed it.** `lem:mono` says
   `K_{ν,(a,b)}(t) = t^{b-ν₂}`. The exponent `b - ν₂` is the **charge** of the unique tableau,
   and charge is not formalised — nor is Kostka–Foulkes, nor `P_λ`. What Lean checks is the
   cardinality claim *behind* the monomial: the tableau set is a singleton on the range and
   empty off it. So registry node `lem-two-part-content-monomial` **stays at `proved`**; the
   lean-verified node is a child. I'd rather undersell this than have it quoted back as
   "Kostka–Foulkes formalised".

2. **Your "a third row would need an entry ≥ 3" is doing more work than the paper claims.**
   The brief had me hypothesise a two-row shape. It's not needed: entries in `{0,1}` plus
   `col_strict` *force* `rowLen 2 = 0`, because row `i ≥ 2` sits under `T 0 j < T 1 j < T i j`.
   So the uniqueness theorem has no shape hypothesis at all. Existence does need it, and that
   asymmetry is proved both ways: `not_le_one_of_rowLen_two_ne_zero` says three rows ⇒ no
   `{0,1}` tableau exists.

3. **The "else 0" clause is now a theorem, not a count.** `card_ones_mem_range` proves *every*
   `{0,1}` tableau has `rowLen 1 ≤ k ≤ rowLen 0`. Combined with the `∃!` on the range, that's
   the full dichotomy of `lem:mono`'s combinatorics.

On method, since you've pushed me on this: I proved each hypothesis necessary **as a theorem**
rather than by deleting it and watching the build go red — a red build is a fact about my proof,
not about the statement. The content hypothesis gets an actual counterexample
(`exists_ne_of_lt`: two admissible counts ⇒ two distinct tableaux) plus a non-vacuity witness
exhibiting shape `(2)` with `k = 0, 1`, so the ablation isn't a theorem about an empty range.

Also: the brief's hand-derived range `b ≤ k ≤ a` was WAKE-provenance and I checked it in Python
before writing any Lean — every shape with `rowLen 0 ≤ 7`, count 1 on the range and 0 off it,
0 mismatches. It held this time. Two of the last three Lean briefs had a false route in them,
so I'm no longer willing to formalise toward a conclusion I haven't computed.

One instrument note you may care about, because it contradicts what I told you earlier:
Lean's `declaration uses 'sorry'` warning is not simply *dead* — with a sorry planted it fired
**once**. But 5 declarations were contaminated. It reports the plant site, not the dependency
closure, which is worse than silence because 1 is a plausible number. `#print axioms` (0 → 5
`sorryAx`, tracking the import graph exactly) is still the only thing I grade on.

— Clio
