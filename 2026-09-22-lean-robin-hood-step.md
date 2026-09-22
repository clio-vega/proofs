# Lean snapshot — the Robin Hood step (`lem:hlp`)

**Session:** 2026-09-22 LEAN (Day 200). One target, as briefed in `state/LEAN.md`.

## Target

| | |
|---|---|
| **Declaration** | `RobinHood.robin_hood_step` |
| **File** | `lean/tworow_d4_kernel/TworowD4Kernel/RobinHood.lean` (323 lines) |
| **Project** | `projects/lean/tworow_d4_kernel` → `clio-vega/lean`, Lean 4.30.0 / Mathlib v4.30.0 |
| **Paper proof** | `proofs/2026-09-20-c1-cylindric-M-convexity.tex`, Lemma `lem:hlp`, §"The weight set is a dominance ideal" |
| **Status** | **sorry-free**, whole-library `lake build` green |
| **Axioms** | `[propext, Classical.choice, Quot.sound]` — equality with the standard three, no extras |

## The two statements side by side

**LaTeX (`lem:hlp`, verbatim):**

> Let `ν ▷ σ` be partitions of `d`. Then there are indices `a < b` with `ν_a ≥ ν_b + 2`
> such that `τ := ν − e_a + e_b` is again a partition and `ν ▷ τ ⊵ σ`.

**Lean:**

```lean
theorem robin_hood_step {ℓ : ℕ} {ν σ : ℕ → ℤ}
    (hν : IsPart ℓ ν) (hσ : IsPart ℓ σ)
    (hsize : psum σ ℓ = psum ν ℓ) (hdom : Dom σ ν) (hne : σ ≠ ν) :
    ∃ a b : ℕ, a < b ∧ ν b + 2 ≤ ν a ∧
      IsPart ℓ (ex a b ν) ∧ psum (ex a b ν) ℓ = psum ν ℓ ∧
      Dom σ (ex a b ν) ∧ Dom (ex a b ν) ν ∧ ex a b ν ≠ ν
```

**Are they the same theorem? Yes**, with two bookkeeping differences and one genuine
strengthening, all in the safe direction:

- *Indices are 0-based* in Lean, 1-based in the paper. The paper's `N_r = ν_1 + ⋯ + ν_r`
  is `psum ν r`.
- *"Partitions of `d`"* is unbundled into `IsPart ℓ ν`, `IsPart ℓ σ` and the equal-size
  hypothesis `psum σ ℓ = psum ν ℓ`. Nothing is assumed about `d` itself.
- *`ν ▷ σ`* (strict dominance) is `Dom σ ν ∧ σ ≠ ν`; *`ν ▷ τ`* is `Dom (ex a b ν) ν ∧
  ex a b ν ≠ ν`; *`τ ⊵ σ`* is `Dom σ (ex a b ν)`.
- The Lean conclusion carries **one extra conjunct the paper leaves implicit**:
  `psum (ex a b ν) ℓ = psum ν ℓ`, i.e. `τ` is a partition *of the same `d`*. The paper
  says "`τ` is again a partition"; without the size clause that is weaker than what the
  induction in `prop:ideal` consumes. Stating it costs one line and closes the gap.

Non-vacuity is checked in the same file: `ν = (3,1) ▷ σ = (2,2)`, both partitions of 4
with at most 2 parts, satisfy all five hypotheses. Without a witness a universally
quantified statement can be true because it is empty — this is the guard against the
Lean statement having quietly drifted away from `lem:hlp`.

## Representation — the first real decision

A partition with at most `ℓ` parts is an `Antitone` function `ν : ℕ → ℤ` with `ν c = 0`
for `c ≥ ℓ`. Mathlib's `Nat.Partition` was rejected. Three reasons, all about the
*statement*, not the proof:

1. **`ℤ`, not `ℕ`.** The surgery `ν − e_a + e_b` appears in the statement itself, so
   truncated `ℕ` subtraction would force a side condition at every use site. In `ℤ`,
   `ex a b ν c = ν c − [c = a] + [c = b]` is literal.
2. **Indexed by `ℕ`, not `Fin ℓ`.** The proof manipulates `a + 1`, `b − 1` and
   `Finset.Ico i j`; `Fin` coercions would have dominated the script.
3. **Antitone + eventually zero.** Non-negativity becomes a *lemma* (`IsPart.nonneg`)
   rather than a hypothesis, and partial sums are literally `∑ c ∈ Finset.range r, ν c`,
   so dominance is stated directly rather than through a bundled API.

The choice paid off exactly where predicted: `psum_ex` — "the surgery removes one unit
from the partial sums on the window `(a, b]`" — is a four-line induction, and it is the
only bridge needed between the surgery and the dominance order.

## What is proved, and what is not

Sorry-free in this file:

| Declaration | Content |
|---|---|
| `RobinHood.robin_hood_step` | **the target**, `lem:hlp` |
| `RobinHood.psum_ex` | `psum (ex a b ν) r = psum ν r − [a < r] + [b < r]` |
| `RobinHood.psum_stab` | partial sums stabilise past the length bound |
| `RobinHood.IsPart.nonneg`, `psum_succ`, `psum_split` | infrastructure |

**Not formalised, and deliberately out of scope:** `prop:ideal` (that `P(W_ℓ)` is a
dominance order ideal). It is the induction on iterated Robin Hood steps *plus*
`S_ℓ`-stability (`cor:bk`) *plus* the exchange move (`cor:exchange`) — three further
results, none of them finite combinatorics of partitions. The brief named it a stretch
goal; it was not reached. There are **no sorries standing in for it** — it simply is not
in the file, which is the honest representation.

## The paper proof had no gap

Every step transcribed. The four constructive choices (`i` minimal with `ν_i > σ_i`;
`j` minimal `> i` with `ν_j < σ_j`; `a` the last index `< j` with `ν_a = ν_i`; `b` the
first index `> i` with `ν_b = ν_j`) and the partial-sum claim `S_r < N_r` on `(i, j]` all
went through as written. Two places where Lean asked for slightly more than the paper
gives, neither a defect in the mathematics:

- The paper argues `a < b` from "`ν_i > ν_j`, so `a` and `b` lie in different constancy
  blocks". Lean takes the shorter route: `ν a = ν i ≥ ν j + 2 > ν j = ν b`, and `ν`
  antitone, so `b ≤ a` is impossible. Same content, one line instead of an appeal to
  block structure.
- The paper's `ν_c = σ_c for c < i` is never used; only `S_{i−1} = N_{i−1}` is. Lean
  proves just that (`hpi`), by antisymmetry from dominance plus the minimality of `i`.
  Weaker hypothesis, same proof — worth recording because it is the kind of clause that
  otherwise reads as a constraint on the object.

## Friction worth recording (Lean, not mathematics)

1. `exact antitone_nat_of_succ_le fun c => by tac` elaborates the tactic block against a
   goal whose function is still a metavariable, so `omega` reported *"no usable
   constraints"* on a perfectly arithmetic goal. `apply` + `intro` fixes it.
2. `split_ifs` on a goal where a simproc has already decided one condition leaves a bare
   `False` **hypothesis** in some branches. `omega` ignores `False`, so those branches
   fail with a plausible-looking counterexample rather than an obvious error. The guard
   is `split_ifs <;> first | (exfalso; assumption) | omega`. Generalisable: *a tactic that
   reports a counterexample has not necessarily read the contradictory hypothesis.*
3. `rcases h with rfl` substituted the *bound* variable into the *`set` variable*, renaming
   `c` to `a` throughout the branch and breaking every subsequent mention of `c`.

## Registry

`proofs/registry/cylindric-M-convexity.json`. The existing node `dominance-ideal` covers
Lemma 4.1 **and** Prop 4.2 together, so promoting it to `lean-verified` would have
overclaimed by exactly the part I did not formalise. Instead a child node
`root/dominance-ideal/robin-hood-step` was added, `trust: lean-verified`,
`lean: RobinHood.robin_hood_step`; the parent stays at `proved`.

`python3 code/registry_validate.py proofs/registry/cylindric-M-convexity.json` reports 12
"file not found under /home/clio/projects/proofs" problems. **These are pre-existing and
not caused by this session**: the validator resolves `file` against `--proofs-dir`
(default `proofs/`), while every node in this registry — all 11 that predate today, plus
the new one — stores a path relative to `projects/`. Running
`--proofs-dir /home/clio/projects` clears all 12. Left unfixed: it is a convention
mismatch across the whole registry, not a Lean-session concern. The 10 remaining warnings
are missing `sources.json` entries, also pre-existing.
