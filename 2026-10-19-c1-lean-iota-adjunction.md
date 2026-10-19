# LEAN 2026-10-19 c1 — the ι-adjunction scalar

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**File:** `TworowD4Kernel/IotaAdjunction.lean` (namespace `IotaAdjunction`)
**Paper proof:** `proofs/2026-10-18-c2-gauge-theorem-and-q-divisors.tex`,
`lem:I2` (L409–L437) and Step 5 of `thm:gauge` (L444–L536).
**Toolchain:** `leanprover/lean4:v4.30.0`, Mathlib `v4.30.0`.

## Target declaration

```lean
theorem iotaAdjunctionScalar (q : K) (y : Fin d → K) (τ : Equiv.Perm (Fin d))
    (hy : ∀ j, y (τ j) = (y j)⁻¹) (hQ : Q2sq q y ≠ 0)
    (hd : ∀ j, 1 - q * (y j) ^ 2 ≠ 0) :
    (1 - q) ^ d * Q2sq q y / Q1 q y = ∏ j, phi q (y j)
```
with `phi q v = (1-q)/(1-q*v^2)`, `Q1 q y = ∏ j, ∏ k, (1 - q * y k / y j)` (all ordered
pairs, diagonal included — paper L388) and `Q2sq q y = ∏_{offDiag} (1 - q * y p.1 * y p.2)`
(the paper's own ordered form of `Q₂²`, L424).

## Status: SORRY-FREE

`#print axioms` — all seven theorems, **standard three only**:

```
'IotaAdjunction.Q1_eq_Gfull'                       [propext, Classical.choice, Quot.sound]
'IotaAdjunction.Gfull_eq'                          [propext, Classical.choice, Quot.sound]
'IotaAdjunction.lem_I2'                            [propext, Classical.choice, Quot.sound]
'IotaAdjunction.lem_I2_div'                        [propext, Classical.choice, Quot.sound]
'IotaAdjunction.iotaAdjunctionScalar'              [propext, Classical.choice, Quot.sound]
'IotaAdjunction.Q1_ne_Gfull_without_iotaClosure'   [propext, Classical.choice, Quot.sound]
'IotaAdjunction.lem_I2_fails_without_iotaClosure'  [propext, Classical.choice, Quot.sound]
```

**sorryAx plant/spare control, both arms reported:**
* planted arm — `theorem plantedSorry : (1:Nat) = 2 := by sorry` → `depends on axioms: [sorryAx]` ✓ fires
* clean arm — all 7 declarations above → `sorryAx` **absent** ✓ spared

**Declarations enumerated from the ENVIRONMENT** (`run_cmd` over `env.constants`), not a
source regex: **18**. Four defs (`phi`, `Q1`, `Q2sq`, `Gfull`), seven theorems, and
**seven auto-generated `iotaAdjunctionScalar._simp_1_1 … _simp_1_7`** emitted by
`field_simp`. A source regex would have reported 11 and missed all seven — the 10-18
`.congr_simp` lesson, reproduced exactly.

## What is lean-verified vs merely cited

**lean-verified:** `lem:I2` and Step 5 of `thm:gauge` — the reindexing identity and the
production of `φ`'s one-variable shape.

**NOT formalised, and this is the honest boundary:** Steps 1–4 of `thm:gauge` (the
`W^C`/`E^C` bookkeeping, the six lines of the `T`-operator, head selections, `cor:heads`,
`lem:head`, `lem:I1`). Nothing in the Lean file mentions the `T`-operator or `r`. So
`thm:gauge` as a whole is **not** lean-verified; its algebraic core is.

## Task 2 — the lesions

The paper's `thm:gauge` carries four hypotheses. The identity spends **exactly one**.

| hypothesis | verdict | evidence |
|---|---|---|
| `y ∘ τ = y⁻¹` (ι-closure) | **IRREMOVABLE** | explicit witness, below |
| `τ` involutive | **REMOVABLE — free strengthening** | never stated; `Q1_eq_Gfull` uses `τ` only as a bijection |
| values in `Fˣ` (`y j ≠ 0`) | **REMOVABLE — free strengthening** | `div_eq_mul_inv`+`inv_inv` is unconditional in a field; at `y i = 0` both sides read `1 - 0` |
| **pairwise distinct** | **NOT SPENT — RELOCATED, not confirmed** | see below |
| `|R| ≥ r` | not expressible | `r` does not occur in the identity |

### The irremovable one, with a witness in a non-empty domain

```lean
theorem Q1_ne_Gfull_without_iotaClosure :
    Q1 (1 : ℚ) (fun _ : Fin 1 => 2) ≠ Gfull (1 : ℚ) (fun _ : Fin 1 => 2)
```
`d = 1`, `y = (2)`, `q = 1` over `ℚ`: `Q₁ = 1 - 1·2/2 = 0` but the full ordered product is
`1 - 1·2·2 = -3`. Domain non-empty (`q = 1 ≠ 0`, `y 0 = 2 ≠ 0`) — the 10-17 `𝔽₂` trap
checked and avoided. Same witness kills `lem:I2` (`lem_I2_fails_without_iotaClosure`).

### `pairwise distinct`: the pre-registration is RELOCATED, not confirmed

Pre-registered: *"it is IRREMOVABLE, and the witness is two coincident deleted values at
the smallest `r ≥ 1`."* What Lean says: **no declaration in the file mentions injectivity
of `y`**, and the identity is proved without it. The hypothesis is not an ingredient of
the reindexing identity at all. In the paper it is attached **upstream**, to
`lem:head`/`cor:heads` (L297, L328), where it makes the line-2a/3a denominator
`(v_p - v_{p'})` nonzero — i.e. it is a **well-definedness hypothesis for the
`T`-operator**. `rem:distinct` already says this in words: the conclusion fails there
*"not because the mathematics is different but because the expression is undefined."*

So the measured `12/12` green at `r = 0` and `0/14` at `r ≥ 1` is a fact about `T`, which
is not formalised here. The retraction was **necessary** for `thm:gauge` as stated — but
the mechanism sits one lemma upstream of Step 5, and Lean cannot see it from here.
**Verdict: the hypothesis belongs to `lem:head`, not to the identity.** This is the same
shape as 10-17's `CollapseMultiplier` finding — *a genericity hypothesis belongs where the
distinction is made, not where the object is described.*

This is also the removable/irremovable pair the brief asked for, and it came out the way
10-18's `a ≠ b` did: **opposite verdicts in adjacent statements**, `hy` irremovable and
`τ`-involutivity removable, neither visible from the other.

## Task 3 — the whole-project build

10-18 c2 recorded `lake build` terminated at `[3463/3471]` with `TworowD4Kernel.olean`
absent. Confirmed at session start: the olean was still missing.

**Diagnosed cause: build-lock contention, and I reproduced it.** Running `lake build`
concurrently with `lake env lean` on 4 cores stalled the full build at `[3468/3471]` for
11 minutes with no error and no completion — the exact 10-18 signature. After killing all
`lake`/`lean` processes and serialising, `lake build TworowD4Kernel.IotaAdjunction`
reported **`Build completed successfully (938 jobs)`, exit 0, in 2.6 s**. The whole-project
build was relaunched alone; its result is recorded in `/tmp/lake-full-2.log`.

**Artifact grading, by SET not by total** (the log cannot be used: lake only prints
modules that emit *warnings*, so absence in the log means "no warnings", not "not built" —
a detector with no positive control, caught and discarded mid-session):

* sources `TworowD4Kernel/*.lean`: **74**; oleans: **78**
* sources with no olean: **{IotaAdjunction}** (mine, now built)
* oleans with no source: **{HalfspaceOrderLesion, ScratchAxioms, ScratchNu, ZZAxiomCheck,
  ZZProbeAutoImplicit}** — 5 stale leftovers from earlier lesion/probe sessions
* `74 - 1 + 5 = 78` ✓ the totals reconcile *only* through the set difference

## Instrument notes

* ⚠ **`pgrep -f "lake|lean"` matched my own session's prompt text** (it contains the string
  "lake build"), reporting two phantom competing processes. Replaced with
  `ps -eo pid,comm | awk '$2=="lake"||$2=="lean"'` — keyed on the *executable*, which
  cannot describe itself — and gave it a positive control (it saw a live `sleep`).
* ⚠ **`Mathlib.Algebra.BigOperators.Ring` does not exist** in Mathlib v4.30.0, and
  `div_inv_eq` is not a lemma name. Both are *predictions about the library, not an
  inventory* — the same class of error as the SageMath line in CLAUDE.md.
* Of the lesion/edit failures this session, **3 of 5 were script-level**
  (`Unknown identifier` from appending past `end IotaAdjunction`; `unknown tactic` from a
  missing `import Mathlib.Tactic.Ring`; `unknown constant Finset.prod_product'` from a
  missing import). Those grade the script, not the statement. Only the two
  `norm_num` witnesses and the `hy` deletion are mathematical.

## Question numbers

**No new question numbers allocated.** Allocator read first-hand at **530**
(`memory/questions/NEXT-FREE.md`) and left at 530 — nothing opened, so nothing bumped,
and no `Q<n>` appears in any identifier in the Lean file (the 10-18 `Q515Plant` fault).
