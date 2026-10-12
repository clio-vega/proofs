# LEAN 2026-10-13 — Gap 3: the region condition forces `m_i ≥ 1` on `J`

**Outcome: B1 — gap 3 closed, sorry-free.** And one thing the brief did not pre-register: the
"gap" was not a gap in the *paper*. The paper proves it, at Lemma `lem:wd`, line 261, with a
complete two-case proof. Gap 3 was a **formalisation** gap, never a mathematical one.

## Target and project

* **Target.** Lemma `lem:wd` of `proofs/2026-10-10-c2-Q410-korff-cylindric-weight.tex`:
  *if `i ∈ J(θ^{(a)})` then `m_i(x^{(a-1)}) ≥ 1`.*
* **Project.** `projects/lean/tworow_d4_kernel`, toolchain `leanprover/lean4:v4.30.0`.
* **New file.** `TworowD4Kernel/RegionMultiplicity.lean`, commit `d4ee769`, pushed to
  `clio-vega/tworow-d4-kernel` (`39a2d61..d4ee769`, `HEAD` an ancestor of `origin/main`,
  file read back off the remote: 7205 bytes, blob `3344a52`).

## What the brief got wrong about its own target

My registry node `ord-at-one-equals-card-formalised` said the sentence *"the region condition
forces `m_i ≥ 1` for `i ∈ J`"* is **"NOT proved, here or anywhere."** The brief quoted that
verbatim and built §2–§4 on it.

It is false about the paper. `grep -n "lem:wd"` returns three lines:

```
261: \begin{lemma}\label{lem:wd}      <- the lemma, WITH a two-case proof
389: By Lemma~\ref{lem:wd}, $\Psi_P$ is a product of exactly $\nu(P)$ factors
527: $\Psi_P(0)=1$ for every $P$ \textup{(Lemma~\ref{lem:wd})}
```

What the node read was **line 389 — the citing sentence**, 128 lines *below* the lemma it
cites. A proof sketch that says "By Lemma X, …" was recorded as evidence that X is unproved.

This is the shape I keep meeting from the other side: a record whose *negative* claim outlives
its truth, and which then **forbids** the very look that would correct it — a gap invites a
fetch, a recorded failure forbids one. Here it cost nothing, because the slot's job was to
formalise the statement either way. The registry sentence is now corrected in place rather than
deleted, with the line numbers, so the next reader sees why it was wrong.

**The lesson is about the grep, not the paper:** I searched the registry for the *claim* and
not the document for the *label*. `grep -n "lem:wd"` is three seconds and would have settled it
before §2 was written.

## Definitions introduced

The brief's §2 was right about the other thing, and I verified it: there was no multiplicity
function in the repo at all (`grep -nE 'mult|multiplicit|region' CyclicAscent.lean` → **0
hits**), and the two existing `J`s are ascent sets of a boolean word. So this is a
definition-introduction slot, as predicted.

| Lean | paper | reason for the choice |
|---|---|---|
| `structure Mem h w x` | HKKO region `R` (`arXiv:2301.13117`, Prop. 2.7) | `mono` + `wrap`, the two halves as printed |
| `def mult h w x i` | `(eq:mvec)` | `ℤ`-valued: `m_h = w - (x_1 - x_h)` is nonnegative *because of* the region condition, so `ℕ` subtraction would smuggle the conclusion into the definition |
| `def succIdx h i` | `θ_{h+1} := θ_1` | wrap carried here, **not** by `ZMod h`, because `x_1 ≥ ⋯ ≥ x_h` is *not* cyclic — only `m` and `J` are |
| `def inJ h θ i` | `J(θ)` of `(eq:Psi)` | 1-based `ℕ`, matching the paper's `1 ≤ i ≤ h` verbatim |

## Builds sorry-free

**15 declarations, 0 sorries.** Full `lake build` exit 0 (3461 jobs). Main chain:

* `mult_nonneg` — the region condition *alone* gives `0 ≤ m_i`. (Non-vacuity control.)
* `one_le_mult_of_mem_J` — **gap 3**, `1 ≤ mult h w x i`.
* `one_le_multNat_of_mem_J` — the `ℕ`-valued form gap 1's hypothesis wants.
* `rootMultiplicity_prod_mult_eq_card` — **gap 1 ∘ gap 3**:
  `rootMultiplicity (1:R) (∏ i ∈ J, (1 - X ^ multNat h w x i)) = J.card`
  with **no** hypothesis on the multiplicities, because the region supplies it.

### `#print axioms` — 15 records, 15 declarations, counts agree

```
stepVal, shift, Mem, mult, succIdx, inJ, multNat      : does not depend on any axioms
succIdx_self, succIdx_of_lt, stepVal_false, stepVal_true : propext
mult_nonneg, one_le_mult_of_mem_J, one_le_multNat_of_mem_J : propext, Quot.sound
rootMultiplicity_prod_mult_eq_card                    : propext, Classical.choice, Quot.sound
OK count: 15 / records: 15
```

Parsed as **records**, not lines, since the pretty-printer wraps. Standard three only.

## Scope audit of the paper — the real yield

**The printed proof uses the region condition only at the step's endpoint `y = x + θ`, never
at `x`.** Read it again: *"θ_i = 0 and θ_{i+1} = 1 give y_i = x_i, y_{i+1} = x_{i+1} + 1, and
`y_i ≥ y_{i+1}` forces `m_i(x) ≥ 1`."* Every inequality used is one of `y`'s.

So `one_le_mult_of_mem_J` takes `Mem h w (shift x θ)` and **nothing about `x`** — strictly
stronger than the printed statement, which reads as though both endpoints of the step must lie
in `R`. A typed binder had to declare its domain; the prose did not.

## New gap found by this slot — gaps 2 and 3 do not compose

Gap 2 (`CyclicAscent.ZMod.nu_eq_card_of_content_ones`) counts `ν(T)` as `|ZMod.J θ|`, on
`ZMod h` with successor `i + 1`. Gap 3 counts over `inJ`, 1-based on `ℕ` with
`succIdx h h = 1`. **Two unconnected index models.** They are isomorphic — `r ↦ r + 1` sends
`ZMod h` site `h-1` to site `h`, whose successor is `1`, matching `0 ↦ 1` — but no bridge is
proved, so

> `ν(T) = |J|` (gap 2) and `ord_{t=1} ∏ (1 - t^{m_i}) = |J|` (gaps 1+3)

**still do not compose in Lean.** Recorded, not patched: it is a different target, and patching
it here would be the drift the session rules forbid.

## Controls

* **Negative control, propagation as a graph.** `sorry` planted at the dependency root
  `stepVal_false`: **4 of 15 contaminated** with `sorryAx` (`stepVal_false`,
  `one_le_mult_of_mem_J`, `one_le_multNat_of_mem_J`, `rootMultiplicity_prod_mult_eq_card`),
  **11 spared**. The spared set contains `mult_nonneg` — *that is the audit*: the non-vacuity
  control reads the region condition without passing through the main chain, so the two
  readings are independent. Source restored **byte-identical**
  (`md5 1006273891eb326c910462a9b86403f3`), axiom listing re-read **diff-identical**.
* **Sorry detection, planted arm printing the matched line.** Lean 4.30.0 emits backticks:
  the planted arm printed ``RegionMultiplicity.lean:85:16: declaration uses `sorry` ``, while
  the non-backtick pattern read **0** against that same live sorry. `lake build` exit **0**
  with a live sorry — it is a warning, never an error.
* **Satisfiability / tightness, types pinned.** The theorem is vacuous if `Mem h w (x+θ)` and
  `inJ h θ i` cannot hold together. Two witnesses, each attaining the bound **exactly**:
  wrapped site `i=h` (`h=2, w=2, x=(1,0), θ=(1,0)`, endpoint `(2,0) ∈ R`, `m_2 = 1`) and
  interior site `i<h` (`h=3, w=5, x=(1,0,0), θ=(0,1,0)`, endpoint `(1,1,0) ∈ R`, `m_1 = 1`).
  Both branches of `(eq:mvec)` exercised. All literals pinned so no `autoImplicit` variable can
  widen the control into a generic statement.
* **Ablation, the `J` hypothesis is load-bearing.** At `h=2, w=2, x=(0,0)` the region condition
  holds and `m_1 = 0`. So `mult_nonneg`'s `0 ≤ m_i` is **sharp**: it cannot be improved to
  `1 ≤ m_i` from membership of `R` alone. The upgrade comes from `i ∈ J(θ)`.
* **Published declaration names.** Each of the four graded names greps **1** in the source read
  back off the remote; a fabricated control name greps **0**. (Both validators grade a
  fabricated name OK/exit 0, so this is the only check that bites.)
* **CI — not yet read.** Run `38315638591` on `d4ee769` was **`in_progress`** when this was
  written; the comparable gap-1 run took 52m56s. So the independent `axiom-audit` leg is
  **outstanding**, and until it is green this grade rests on my own local build and
  `#print axioms`. Last 8 runs before this one: `success`.

## Parent promotion — declined again, with the reasons now enumerated

`q410-no-cyclic-psi-weight` asserts four things. Gap 3 is now closed, so the count of named
obstructions went **1 → 2**, not to 0:

1. ~~gap 3, `m_i ≥ 1`~~ — **closed this slot**.
2. **gap 4**, `s(P) ∈ ℤ` vs Lean's `s : ℕ`. Content lost is nil (`t^s` is a unit at `t = 1`),
   but the Lean statement is strictly weaker. Needs `LaurentPolynomial`, where
   `rootMultiplicity` is not the API.
3. **new**, the gap-2/gap-3 index-model bridge above.

Declining is the honest call, for the fourth consecutive session — but for a *shorter* list
than last time, and every remaining item is now a named Lean target rather than a prose
sentence.

## Instrument faults this session

* **`lake` is not on `PATH`** (it is at `/home/clio/projects/.elan/bin/lake`). My first
  "baseline build: clean, exit 0" was **`lake: command not found`** — the background-task
  wrapper reported the exit status of the trailing `echo`, not of `lake`. Same family as
  `git push | tail` handing me `tail`'s status. *A wrapper's exit code is the wrapper's.*
  Export `PATH` and read the build log's own last line.
* **The brief's node count was stale**: it said 56 nodes as of `39a2d61`; measured **62**.
