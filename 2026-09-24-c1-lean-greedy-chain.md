# LEAN snapshot — 2026-09-24 c1 — the greedy chain

**Project:** `/home/clio/projects/lean/tworow_d4_kernel/`, module `TworowD4Kernel.GreedyChain`,
Lean 4 / Mathlib v4.30.0.
**Repo:** `github.com/clio-vega/tworow-d4-kernel`, commit **`1d35afb`**.
**Paper proof:** `proofs/2026-09-20-c1-cylindric-M-convexity.tex`, §"The greedy chain and
the dominance maximum" (`def:greedy`, `eq:greedy`, `lem:greedy`, `prop:noell`).

## Outcome

`lake build` green on the whole project (2989 jobs). **0 sorries** — `grep sorry` finds one
hit, in the module docstring saying there are none. No `native_decide`, no local `axiom`.

All three milestones of the brief landed, plus a non-vacuity witness that was not asked for.

## What builds, sorry-free

| Declaration | Paper statement | Axioms |
|---|---|---|
| `greedy` | `def:greedy` / `eq:greedy` | none |
| `greedy_isCylindric` | `lem:greedy`(1), first sentence — **A1** | `propext, Quot.sound` |
| `greedy_sub_lam` | `gᵗ ⊆ λ` | `propext` |
| `greedy_le_succ` | `gᵗ ⊆ gᵗ⁺¹` — **C1** | `propext, Quot.sound` |
| `greedy_hstrip` | `lem:greedy`(1) — **A2** | `propext, Quot.sound` |
| `chain_mono` | chains are weakly increasing | `propext, Quot.sound` |
| `chain_sub_lam` | `xᵗ ⊆ xˡ = λ` | `propext, Quot.sound` |
| **`greedy_dominates`** | **`lem:greedy`(3) — B** | **`propext, Quot.sound`** |
| `u_le_u_greedy` | `lem:greedy`(3), weight form | `propext, Classical.choice, Quot.sound` |
| `greedy_fix` | `prop:noell` — **C2** | `propext, Quot.sound` |
| `greedy_eq_lam_of_le` | `gˢ = λ` for `s ≥ ℓ₀` | `propext, Quot.sound` |
| `incr_eq_zero_of_ge` | `γₜ = 0` for `t > ℓ₀` — **C3** | `propext, Classical.choice, Quot.sound` |
| `nonzeroIncr_indep_of_len` | `prop:noell`, `ℓ`-independence | `propext, Classical.choice, Quot.sound` |
| `witness_mu_cylindric` | non-vacuity | `propext, Quot.sound` |
| `witness_lam_cylindric` | non-vacuity | `propext, Quot.sound` |
| `witness_dominates_strictly` | non-vacuity + strictness | `propext, Classical.choice, Quot.sound` |

Standard three throughout, never more. **The whole mathematical spine — A1, A2, B, C1, C2 —
is choice-free**: `[propext, Quot.sound]` only. `Classical.choice` enters exactly where
`Finset.sum` and `Multiset.filter` do, i.e. in the weight-level restatements, not in the
domination theorem.

## Nothing is assumed

Yesterday's `AffineAdditive.lean` had to *define* the paper's Lemma 2.3 criterion because the
affine symmetric group is not in Mathlib — a transport the type checker cannot falsify, which
is why it owed an external differential check.

Today there is no transport. `IsCylindric` (`per` + `inc`), `Sub`, `HStrip` and `greedy` are
literal transcriptions of `def:shape`, `eq:hstrip` and `eq:greedy`; every one is an explicit
inequality on `ℤ → ℤ`. The only unfalsifiable step left is my own transcription from the
`.tex`, and that is what the Python check below is for.

Two things I want on the record about the statements themselves:

* **`greedy_dominates` needs neither periodicity nor `μ ⊆ λ`.** It uses only the strip
  conditions. The paper bundles the hypotheses of `lem:greedy`(1) and (3) together; in Lean
  they separate, and (3) is the weaker statement. (Cf. the standing lesson that a hypothesis
  of my own theorem is an unexamined assumption — here dropping two of them cost nothing.)
* **`greedy_fix` (C2) is one line and it is the whole of the correction to gap (iii).**
  `min(λᵢ₊₁ − 1, λᵢ) = λᵢ` because `λᵢ < λᵢ₊₁`. Strict increase is doing all the work. The
  09-20 note's clause "λ̂ grows with ℓ" was false, and this is the sentence it missed.

## Non-vacuity

`greedy_dominates` is an implication and would be satisfied by an empty hypothesis set, so
`witness_dominates_strictly` exhibits `n = 2`, `m = 1`, `μ i = 2i`, `λ i = 2i + 1` and the
legal chain `μ, μ, λ`, whose middle term satisfies `x¹ 0 = 0 < 1 = g¹ 0`. So the hypotheses
are satisfiable **and** the conclusion is not secretly an equality.

## Differential check

`proofs/code-greedy/greedy_check.py` re-implements `def:shape`, `eq:hstrip` and `eq:greedy`
independently and enumerates **all** chains by brute force (successors are built straight from
`eq:hstrip` as a box, then filtered for cylindricity — `lem:box` is deliberately *not*
assumed, so the check does not inherit the paper's reasoning).

Parameters `2 ≤ n ≤ 5`, `1 ≤ m ≤ 3`, chain length `≤ 4`:

```
pairs        397
a1           4367
a2           3970
b_chains     5925
b_coords     60225
c_fix        2420
c_indep      867
nonempty     1985
fail         0
```

**0 failures.** `b_coords` is the count of coordinate comparisons behind `greedy_dominates`
across every chain enumerated; `nonempty` checks `prop:noell`'s other half.

A wider run (`2 ≤ n ≤ 7`, chain length `≤ 6`) was launched and had not finished inside the
session; it reports failures the moment it finds one and had reported none. The counts above
are from the run that completed, and only those are claimed.

## What is NOT done, and why

1. **`lem:greedy`(2), `gˡ = λ`.** The paper derives it from (3) applied to any chain, under
   the hypothesis `𝒯_ℓ ≠ ∅`. Formalising it means carrying that nonemptiness hypothesis
   around; it is a half-hour of plumbing, not a difficulty, and it was below A/B/C in the
   brief's order. The Python check verifies it (`nonempty`, 1985 instances).

2. **`prop:max`** — that the greedy weight *is* the dominance maximum of the support. Needs
   `Sₗ`-stability and the Bender–Knuth corollary. Explicitly out of scope per the brief.
   `greedy-maximum` therefore stays at `proved` in the registry; only the new child is
   `lean-verified`.

3. **`sort` versus the multiset of nonzero increments.** `nonzeroIncr_indep_of_len` proves
   that the *multiset of nonzero increments* is `ℓ`-independent. The paper says "sorting
   discards the trailing zeros", which additionally asserts

   > `List.filter (· ≠ 0) (Multiset.sort (· ≥ ·) γ) = Multiset.sort (· ≥ ·) (γ.filter (· ≠ 0))`

   — that passing to the sorted presentation commutes with dropping zeros. That is a genuine
   proposition, not presentation, and it is **not proved**, here or on paper. I am recording
   it as a proposition rather than waving it through: a clause filed as cosmetic is an
   unchecked claim. It is true and it is easy (`List.eq_of_perm_of_sorted`); it is simply not
   done, and the brief capped C3 at 40 minutes.

## Registry

`proofs/registry/cylindric-M-convexity.json`: new child of `greedy-maximum`,

* `id`: `greedy-dominates-coordinatewise`
* `trust`: `lean-verified`
* `lean`: `GreedyChain.greedy_dominates`
* `file`: `lean/tworow_d4_kernel/TworowD4Kernel/GreedyChain.lean`

`greedy-maximum` itself **not** promoted. Backup of the whole registry directory at
`proofs/registry.bak-20260924-c1-lean/`.

`trustcheck ... --files-dir .` → `OK: ... is valid`.

`registry_validate.py` prints 14 `file ... not found under /home/clio/projects/proofs`
problems, one per node, **all spurious and all pre-existing**: it resolves `file` fields
against `projects/proofs/` when they are relative to `projects/`. I checked rather than
assumed — `ls lean/tworow_d4_kernel/TworowD4Kernel/GreedyChain.lean` from `projects/`
succeeds. This is the same root-directory defect already recorded for
`trustcheck --files-dir`, in a second tool; the fix belongs in `registry_validate.py`'s
default root and is not mine to make in a Lean session. It also exits 0 while printing
`14 problem(s)`.

I also dropped a `sources` entry I had first written for the new node
(`"2401.14632 Section 5 ..."`): it is absent from the sources index, and minting a fresh
identifier that resolves to nothing is exactly what the index exists to prevent. The node's
content is transcribed from my own `.tex`, which is already in its `file` field, and the Lean
docstring carries the arXiv citation.

## Reproduce

```
export PATH="/home/clio/projects/.elan/bin:$PATH"
cd /home/clio/projects/lean/tworow_d4_kernel && lake build TworowD4Kernel.GreedyChain
cd /home/clio/projects/proofs/code-greedy && python3 greedy_check.py
```

Note for future sessions: `lake` is at `/home/clio/projects/.elan/bin`, **not** `~/.elan`,
and is not on the default `PATH`. `import Mathlib.Tactic.Omega` does not exist in v4.30.0
(omega is core); `interval_cases` needs an import this file does not carry.
