# LEAN 2026-09-23 c1 — `AffineAdditive`: the length-additivity criterion, transported

**Project:** `lean/tworow_d4_kernel` (Lean 4.30.0 + Mathlib)
**Module:** `TworowD4Kernel/AffineAdditive.lean`, imported by the root aggregator
**Paper:** `proofs/2026-09-23-c1-bracketed-pair-deletion.tex` (today's PROVE output —
a complete proof, not a failure report), on `lem:blocks`/`lem:add` of
`proofs/2026-09-21-c1-affine-stanley-exchange.tex`
**Registry:** new node `…/emptyX-two-factor-extension/additivity-not-hereditary`,
`trust: lean-verified`

## Result: 13 declarations, 0 sorries, standard axioms only

`lake build` green. `grep sorry` finds one hit, in a docstring. No `native_decide`,
no local `axiom`.

| declaration | what it is | axioms |
|---|---|---|
| `firstGapAux_periodic` | `n`-periodicity of the gap search | `[propext, Quot.sound]` |
| `firstGap_periodic` | `firstGap S (j+n) = firstGap S j + n` | `[propext, Quot.sound]` |
| `uS_periodic` | `u_S (j+n) = u_S j + n` (eq. `equiv` of the note) | `[propext, Quot.sound]` |
| `criterion_periodic` | the run condition is translation-invariant | `[propext, Quot.sound]` |
| `additive_zero_zeroone` | `({0},{0,1})` additive at `n=3` | standard three |
| `not_additive_zero_zero` | `({0},{0})` not additive at `n=3` | standard three |
| `additive_zeroone_one` | `({0,1},{1})` additive at `n=3` | standard three |
| `not_additive_one_one` | `({1},{1})` not additive at `n=3` | standard three |
| **`additive_not_hereditary`** | **Milestone A** | standard three |
| **`additive_erase_left_not_stable`** | **the corrected Milestone A** | standard three |
| `deletion_n3` | deletion theorem, all `(S,T,i)` at `n=3` | standard three |
| `deletion_n4` | deletion theorem, all `(S,T,i)` at `n=4` | standard three |
| `deletion_n5` | deletion theorem, all `(S,T,i)` at `n=5` | standard three |

"standard three" = `[propext, Classical.choice, Quot.sound]`, verbatim. Nothing else
appears anywhere. The four periodicity lemmas are *choice-free*.

## What is formalised, and what is assumed

The affine symmetric group is **not** built (not in Mathlib, and not what the theorem
is about). `Additive S T` is *defined* to be the criterion of `lem:add` applied to the
block model of `lem:blocks`:

```lean
def uS (S : Finset (ZMod n)) (j : ℤ) : ℤ :=
  if ((j - 1 : ℤ) : ZMod n) ∈ S then j - 1 else firstGap S j

def Additive (S T : Finset (ZMod n)) : Prop :=
  S.card < n ∧ T.card < n ∧
    ∀ k : Fin n, ((k.val : ℤ) : ZMod n) ∈ T →
      uS S (k.val : ℤ) < uS S (firstGap T ((k.val : ℤ) + 1))
```

`lem:blocks` and `lem:add` are **assumed**, exactly as 2026-09-22 assumed `cor:bk` and
`cor:exchange`. So this is a formalisation of the **combinatorial core**, not of the
group-theoretic statement. Said plainly in the module header.

## The transport is checked, not trusted

A defined-not-derived predicate is only as good as the claim that it *is* the thing.
`proofs/code-q232/lean_transport_check.py` re-implements the Lean definitions verbatim
and compares them against `ℓ(u_S u_T) = |S| + |T|` computed by **Shi's inversion
formula** on the composed affine permutation — a route that never touches the block
model:

> **0 disagreements over all 87 376 pairs, `2 ≤ n ≤ 8`.**

Two further guards:

* `Additive` quantifies over lifts `k ∈ {0,…,n-1}` rather than over all runs in `ℤ`.
  The paper licenses this by `n`-equivariance. Here it is **proved**, not assumed:
  `criterion_periodic`.
* `#eval` shows `u_{{0}} = [1,0,2]` (i.e. `s₀`) and `u_{{0,1}} = [2,0,1]` at `n = 3`,
  each of length `|S|` — the block model behaves.

## Finding: the session brief's Milestone A proves the wrong thing

The brief said the `n=3`, `S={0}`, `T={0,1}`, `S'=T'={0}` witness shows "the bracket
condition in Q232 is not decoration", and that a Lean/Python disagreement about it
would be the finding of the session. **Lean and Python agree** — the witness is real
and `additive_not_hereditary` proves it. The finding is elsewhere: *the inference from
that witness does not go through*, and today's PROVE note says so independently
(§"The bracket hypothesis is not load-bearing").

In that witness **nothing is deleted from `S` at all**. The deletion is `T ↦ T∖{1}`
with `S` held fixed. So it says nothing about whether `i+1 ∈ T` is needed. Indeed the
paper proves the theorem *without* that hypothesis, in a stronger uncoupled form.

I therefore formalised both statements and named them for what they are:

* `additive_not_hereditary` — the brief's statement, true, and labelled in its
  docstring with what it does not show.
* `additive_erase_left_not_stable` — the claim that **is** load-bearing: deleting `i`
  from `S` while leaving `T` alone breaks additivity (`n=3`, `S={0,1}`, `T={1}`,
  `i=0`, giving `({1},{1})`). This is the paper's *second* witness, and it is the one
  that shows the two deletions must be **coupled**.

`deletion_n3/n4/n5` then machine-check the theorem itself — in the strengthened form,
hypothesis `i+1 ∈ T` absent — exhaustively over every `(S,T,i)` at those `n`. `n=5` is
the first size at which the `ms-crystal-comparison` containment is strict, so it is
the first `n` not forced by smaller cases.

## Milestone B: not attempted, and why

The general-`n` proof is **not** formalised, and nothing stands in for it — no `sorry`,
no `axiom`. The paper proof needs `lem:perturb` (deleting `i` from `S` moves `u_S` at
exactly two positions, in opposite directions), which with the fuel-bounded `firstGap`
definition means proving the block-splitting `[p,q+1] ↦ [p,i] ⊔ [i+1,q+1]` from
scratch, plus `lem:below`, plus the three-case integer-interval analysis. That is a
multi-session job, not a 20-minute one. Recorded as the next target rather than sorried.

**Which case does the paper wave at?** None, as far as this session can tell — and I
looked for the predicted "obvious" side condition on `i` being the bottom or top of its
run. It is not there: `lem:incong` (`p ≢ i+1`) is exactly that side condition, and the
paper *does* state and prove it, from `q - p ≤ n-2`. The one place the prose moves fast
is inside Case 1b, the parenthetical "and if `p < m'` there is no such `j` at all,
which is also fine" — a vacuous sub-branch asserted rather than argued. It is correct,
but Lean would make it an explicit case.

## Reproduce

```sh
. /home/clio/projects/lean/ENV.sh
cd /home/clio/projects/lean/tworow_d4_kernel && lake build TworowD4Kernel.AffineAdditive
python3 /home/clio/projects/proofs/code-q232/lean_transport_check.py
```
