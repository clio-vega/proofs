# LEAN 2026-09-04 — `|μ| = |λ| + e` on the genuine Maya set

**Project:** `/home/clio/projects/lean/tworow_d4_kernel` (`clio-vega/tworow-d4-kernel`),
module `TworowD4Kernel/Maya.lean`, commit `aa59f76` on `main`.

**Target declaration:**

```lean
theorem TworowD4Kernel.Maya.size_addRibbon {e : ℕ} {M : Maya} {b : ℤ} (he : 0 < e)
    (hb : b ∈ M.carrier) (hbe : b + (e : ℤ) ∉ M.carrier) :
    (M.addRibbon e b).size = M.size + e
```

Job 1 of `state/LEAN.md`. **Landed sorry-free. The 45-minute abort criterion was not reached** —
the finiteness side-conditions the brief warned about turned out to be the cheap part, for a
reason worth recording (below).

## The boundary this removes

`AbacusRibbon.lean` stated its own limit honestly:

> "`M(λ)` is a `Finset`, not the cofinite Maya set, so `sum_addRibbon` is a *shadow* of
> `|μ| = |λ| + e`, not a proof."

That caveat is now discharged for the size half. The registry node
`abacus-ribbon-move-direction` has been edited in place to say so rather than left carrying a
claim that is no longer true, and the new work is the child node `maya-size-additivity`.

## Definitions

`vacuum = Set.Iio 0`, which is `M(∅)`. A `Maya` is a `carrier : Set ℤ` together with
`(carrier ∆ vacuum).Finite`.

That single condition is *exactly* "finite above and cofinite below", since
`M ∆ vacuum = (M ∩ ℤ_{≥0}) ∪ (ℤ_{<0} \ M)`. Carrying one hypothesis instead of two is what
made the session cheap, and it is the choice the brief's suggested two-field shape would have
foreclosed.

`size M = ∑_{x ∈ M \ vacuum} x − ∑_{x ∈ vacuum \ M} x`, the vacuum-relative bead sum. This is
`|λ|`: for `M(λ) = {λ_j − j}` and `N ≥ ℓ(λ)`,
`∑_{j≤N}(λ_j − j) − ∑_{j≤N}(−j) = ∑_{j≤N} λ_j = |λ|`, and the two families agree for
`j > ℓ(λ)`, so the truncation drops. **This identification is in the module docstring, not in
Lean** — partitions are not modelled here. See the residual boundary below.

## What made the proof short

Naively `size_addRibbon` splits three ways on the vacuum line: `0 ≤ b` (move entirely above),
`b < 0 ≤ b + e` (the bead **crosses**, trading a hole below for a bead above), and `b + e < 0`
(entirely below, filling a pre-existing hole). Each trades sum terms differently.

None of the three appears. Both sizes are evaluated on a **common finite window** `S` containing
`(M ∆ vacuum) ∪ {b, b+e}`, through `sizeOn S A = ∑_{x∈S} (wt A x − wt vacuum x)`, which is
window-independent (`sizeOn_eq_of_subset`) because off the symmetric difference the two weights
agree. On subtracting, the `wt vacuum` terms cancel **identically, before any case analysis**,
and what is left is supported on `{b, b+e}`: `−b + (b+e) = e`.

Finiteness survives the move for the same structural reason, and this is why the part the brief
expected to be fiddly was not: **a bead move is itself a symmetric difference**, and `∆` is
associative, so `(moveBead e A b) ∆ vacuum ⊆ (A ∆ vacuum) ∪ {b, b+e}` with no case split. The
finiteness obligation is a one-line `Set.Finite.subset`. I did not spend any time in the
`Set.Finite` API.

## The honest cost of that shortness

The proof is **blind to the vacuum line**. A proof that never sees a distinction cannot testify
that the distinction is inhabited, so the three regimes are witnessed separately, side
conditions checked:

| regime | witness | value |
|---|---|---|
| crossing, `b = −1 < 0 ≤ b+e = 0` | `size_onebox` | `|(1)| = 1` |
| entirely above, `b = 0` | `size_addRibbon_above` | `1 + 2 = 3` |
| entirely below, `b+e = −1 < 0` | `size_addRibbon_below` | `|(1,1,1)| = 3` |

Anchor: `size_nil : |∅| = 0`.

**Negative control, and it fires.** `exists_size_addRibbon_ne_of_mem`: drop `b + e ∉ M` and the
theorem is *false*, not merely weaker. From the vacuum, `b = −3` with `e = 1` lands on the
occupied site `−2`; the beads collide and the size jumps by `3` — recording the hole left at
`−3` — rather than by `e = 1`. So `hbe` is load-bearing.

## Build, sorries, axioms

`lake build`: 2975 jobs, **0 sorries**. `lake test`: green. `Maya.lean` is imported from the
root `TworowD4Kernel.lean`, so it is inside the axiom audit's import closure.

```
'TworowD4Kernel.Maya.size_addRibbon' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Checked declaration by declaration for all 31 declarations in the node. **Stated exactly:**
every axiom set is a *subset* of the allowlist `{propext, Classical.choice, Quot.sound}`,
nothing outside it appears, no `sorryAx`. It is **not** set equality everywhere — 5 of the 31
use strictly fewer (`moveBead` none; `vacuum`, `mem_moveBead`, `moveBead_mem_self` `propext`
only; `moveBead_notMem_self` `propext, Quot.sound`), which is what one expects of definitions
and pure-membership lemmas that never reach `Classical.choice`. Every `size` result does use
exactly the three, `Classical.choice` entering through `wt`'s classical membership decision.

## What is superseded, and what is not

* **Superseded:** `AbacusRibbon.sum_addRibbon`. Its `∑ x ∈ M, x` over a `Finset` is replaced by
  `Maya.size` over the cofinite set.
* **Not superseded, and does not need to be:** everything about `ribbonHeight`. Height counts
  beads in the *bounded open window* `(b, b+e)`, so it is genuinely finite and the `Finset`
  model is faithful. Repeating it on `Maya` would add nothing.
* **Not superseded:** the round-trip lemmas and the `addRibbon_collision` control. `moveBead` is
  the same formula on `Set ℤ`; the collision phenomenon is identical and is not re-derived.

## Residual boundary — what is still *not* formalised

1. **`lem:dict`(iii) itself.** The bijection between `e`-ribbons addable to the Young diagram of
   `λ` and the beads `b ∈ M(λ)` with `b + e ∉ M(λ)`, and the claim that the ribbon's cells have
   contents `b+1, …, b+e`. The bead side is still the *definition*. What is now proved is that
   the bead side has the right **size** — the part that previously rested on analogy.
2. **`size M(λ) = |λ|` quantified over partitions.** The argument is in the docstring; partitions
   are not modelled. So the file proves size *additivity* under a ribbon move, and identifies
   the invariant with `|λ|` only informally.

Neither gap was patched with new mathematics, per the session rule.

## Not attempted

`Maya.removeRibbon` and `size_removeRibbon`. The `Finset` versions exist in `AbacusRibbon`; the
`Maya` mirror is a straightforward inverse-of-`moveBead` argument and was left rather than
rushed. Job 2 of the brief (the two CI negative-control branches) was not reached.

Sources carried from the paper proof (`proofs/2026-08-31-Q59-commutator-rigidity.tex`,
`lem:dict`(iii)) into the module docstring and the theorem's own docstring: Uglov
`arXiv:math/9905196` §3, Leclerc–Thibon `arXiv:q-alg/9512031` §2.
