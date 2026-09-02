# Lean session 2026-09-02 — the window lemma

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
(GitHub `clio-vega/tworow-d4-kernel`, commit `00c7e6d`)
**Module:** `TworowD4Kernel/WindowLemma.lean`, imported from the root `TworowD4Kernel.lean`

## Target

> In any `N` consecutive integers, each residue class mod `N` occurs exactly once.

This is the arithmetic engine behind `0 ≤ a ≤ ℓ − 1` in the closed-form theorem of
`2026-09-01-Q63-wedge-fock-normalisation.tex` (§`sec:closed`, 913/913). There the `ℓ`
runners of the abacus are the residue classes mod `N` in Uglov's indexation
`k = c + n(d−1) − n ℓ μ` (`math/9905196`), and the locality argument needs a window of
`N` consecutive positions to meet every runner exactly once — no more (which bounds how
many beads can move) and no fewer (which is what makes the bound attained).

## What builds, sorry-free

Both declarations, in namespace `TworowD4Kernel`:

```lean
/-- Indexing form. -/
theorem window_residues_bijective (N : ℕ) [NeZero N] (p : ℤ) :
    Function.Bijective (fun i : Fin N => ((p + i : ℤ) : ZMod N))

/-- Counting form. -/
theorem card_window_eq_one (N : ℕ) [NeZero N] (p : ℤ) (r : ZMod N) :
    ((Finset.Ico p (p + N)).filter (fun k : ℤ => (k : ZMod N) = r)).card = 1
```

plus two `private` helpers, `windowRep` (the canonical representative
`p + (r − p).val` of the class `r` inside the window) and its two properties
`windowRep_mem`, `windowRep_cast`.

Both forms are kept because the downstream argument uses both: the indexing form to
*name* the runner of the `i`-th position, the counting form to *count* beads on a fixed
runner.

**Sorry count: 0.** `lake build TworowD4Kernel` → `Build completed successfully
(2972 jobs)`, exit 0, no `declaration uses 'sorry'` warnings.

## `#print axioms`

```
'TworowD4Kernel.window_residues_bijective' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.card_window_eq_one' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Asserted as **equality** with the standard three, not as a substring check. A substring
check is unsound here: `native_decide` in Lean 4.30.0 emits a per-declaration
`<decl>._native.native_decide.ax_N_M` and never `Lean.ofReduceBool`, so a blocklist keyed
on that name passes a `native_decide` proof silently. Neither proof uses `decide` or
`native_decide` in any case — both are structural.

## Import closure — the precondition, honoured first

`axiom-audit` follows the **import closure of `axiom-audit-root`**, not the namespace. The
`TworowD4Kernel.+` glob in `lakefile.toml` makes `lake build` compile a module; that is
*not* the same thing, and a module reached only by the glob is invisible to the audit.

So the root import line was added **before any mathematics was written**. Recomputed
closure after the change: **19 modules** (was 18), `TworowD4Kernel.WindowLemma` among them.

## Was it already in Mathlib?

Searched `Mathlib v4.30.0` before proving: `Mathlib/Data/ZMod/`, `Mathlib/Order/Interval/Finset/`,
`Mathlib/Data/Int/Interval.lean`, and for `ExistsUnique`/`∃!` over `Ico` with a modular
hypothesis. The **ingredients** are all there —

* `ZMod.intCast_zmod_eq_zero_iff_dvd` (`Data/ZMod/Basic.lean:511`)
* `Int.eq_zero_of_abs_lt_dvd` (`Algebra/Order/Group/Unbundled/Int.lean:83`)
* `ZMod.val_lt`, `ZMod.val_cast_of_lt` (`Data/ZMod/Basic.lean:1003`)
* `Fintype.bijective_iff_injective_and_card`, `ZMod.card`

— but **neither statement exists in this form**. That is a claim about a search, not a
theorem, and it should be re-tested on any Mathlib bump.

## Choices worth recording

**`ZMod N`, not `Int.emod`.** The brief asked for this and it is the right call for a
reason worth writing down: the runner index is an element of a cyclic group of order `N`.
It is *not* a distinguished representative in `[0, N)`. Phrasing the lemma over `% N`
would have made the statement true but would have quietly imported a choice of
representative that the abacus argument does not make, and the downstream use — "the label
`(c', d')` occurs exactly once in the window" — is about the class, not the representative.

**Where the `Int`/`ZMod` bridge actually lives.** Not in a cast lemma; in
`Int.eq_zero_of_abs_lt_dvd`. Congruence gives `N ∣ k − w`; the window gives
`|k − w| < N`; those two together give `k = w`. Both halves of the window lemma are that
one move, dressed differently. The `push_cast`/`ZMod.natCast_val` layer is bookkeeping.

**One elaboration trap, recorded because it cost a round trip.** Writing the filter as
`(fun k => (k : ZMod N) = r)` — the form suggested in the brief — makes Lean read
`(k : ZMod N)` as a *type ascription* on the binder rather than a coercion, and the whole
`Finset.Ico p (p + N)` then fails to unify: `Finset ℤ` against `Finset (ZMod N)`. The
binder must be annotated: `(fun k : ℤ => (k : ZMod N) = r)`.

## What is NOT formalised

The abacus. Uglov's indexation, the identification of the `ℓ` runners with the residue
classes mod `N`, the locality lemma (straightening `k_j → k_j + N` never disturbs a bead
outside the closed window `[k_j, k_j + N]`), and the step from "each label occurs once in
the window" to `a ≤ ℓ − 1` are all still **on paper only**, with 913/913 computational
checks behind them. This module formalises the single arithmetic fact that argument rests
on — which is what turns `a ≤ ℓ − 1` from an observation over `|λ| ≤ 4` into something
with a proof, and no more than that.

## Registry

Node `Q63-window-lemma` added under
`root/Q63-level-ell-telescope/Q63-heisenberg-not-single-bead/Q63-wedge-fock-normalisation/Q63-B-single-component-closed-form/Q63-exponent-range-and-cross-runner`
in `proofs/registry/fock-ribbon-sign-operator.json`, graded `lean-verified`, with both
declaration names in the `lean` field.

`registry_validate.py` reports one remaining problem in this file, **pre-existing and
unrelated**: `Q63-rigidity-dichotomy` claims `proved` while its child
`Q63-formii-empty-at-level-1` is `computed`. Not touched this session.

## Session cost note (for the next Lean session)

`lake build` in this container is **I/O-bound, not compute-bound** — the first invocation
spent ~10 min at ~8 s of user CPU replaying trace hashes over the 6.2 GB of Mathlib
oleans, and timed out twice before completing. Two things that helped and should be reused:

1. Run the full `lake build` in the background once to warm the trace cache, and
2. iterate on the module with `lean` directly, bypassing Lake's lock and trace check:
   `LEAN_PATH=<project>/.lake/build/lib/lean:<each package>/.lake/build/lib/lean lean Module.lean`
   (~40 s per round trip against ~10 min).
