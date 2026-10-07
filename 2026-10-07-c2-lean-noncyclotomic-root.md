# LEAN 2026-10-07 c2 — the real root of `t^b − t^{b−1} + 1`

**Target:** `TworowD4Kernel.exists_root_Ioo`
**Project:** `/home/clio/projects/lean/tworow_d4_kernel/`
**New file:** `TworowD4Kernel/NonCyclotomicRoot.lean` (added to the root aggregator
`TworowD4Kernel.lean`; `TypedBlocks.lean` untouched)
**Paper proof:** Theorem D, `proofs/2026-10-07-two-part-green-polynomials.tex`

```lean
theorem exists_root_Ioo (b : ℕ) (hb : 3 ≤ b) (hodd : Odd b) :
    ∃ t ∈ Set.Ioo (-1 : ℝ) 0, t ^ b - t ^ (b - 1) + 1 = 0
```

**Sorry-free. 8 declarations, 0 sorries.** The must-have landed; the stretch goal
(non-cyclotomicity) was not started — see *What is not formalised*.

---

## Pre-flight guard, run on the statement shape

Ran as mathematics, not as a name:

- `Set.Ioo` — **zero** occurrences in the library outside `Finset.Ioo` (a different object: the
  abacus window code in `AbacusRibbon`/`AdamsDilation`).
- `intermediate_value` — **zero**. `Continuous`/`ContinuousOn` — **zero** in any
  real-analytic sense; the only three `Mathlib.Analysis` imports in 858 declarations are
  `Matrix.Spectrum` (×2) and `Convex.Hull`.
- `IsRoot` — three hits, all finite-field or `Polynomial.reflect` (`Fp2Irreducible`,
  `SelfReciprocal`), none about real roots.
- `Odd b` as a hypothesis — hits in `B0modKernel`, `ThreeRowC4Boundary`, `GaussianUnitSum`;
  read each, all 2-adic or parity-of-a-sum, none about sign at `−1`.

Neighbours whose docstrings I read rather than trusting their names: `PadicNoRoot.lean` is
*p*-adic non-existence for `X²+X+1`; `SelfReciprocal.lean` is multiplicity parity at `t = −1`
via `reflect`. Neither states this. **The slot was genuinely empty.**

Searched Mathlib for what the proof *uses*, not what the statement *says* — and the thing the
statement's noun would have found (`Polynomial.roots`) is not what the proof needs:
`intermediate_value_Ioo (hab : a ≤ b) (hf : ContinuousOn f (Icc a b)) : Ioo (f a) (f b) ⊆ f '' Ioo a b`
at `Mathlib/Topology/Order/IntermediateValue.lean:617` is an exact fit.

## What builds

| declaration | content |
|---|---|
| `D` | `D b t = t^b − t^(b−1) + 1` as a bare `ℝ → ℝ` |
| `continuous_D` | `Continuous (D b)` |
| `D_eval_neg_one` | `Odd b → D b (−1) = −1` — **the only use of the parity hypothesis** |
| `D_eval_zero` | `2 ≤ b → D b 0 = 1` |
| **`exists_root_Ioo`** | **the target** |
| `D_four_eval_neg_one` | `D 4 (−1) = 3` |
| `D_four_pos_of_mem_Icc` | `t ∈ Icc (−1) 0 → 0 < D 4 t` |
| `not_exists_root_Ioo_four` | the ablation *as a theorem* — see below |

Lean mirrors the paper line for line: `Odd.neg_one_pow` and `Nat.Odd.sub_odd` give
`(−1)^b = −1`, `(−1)^(b−1) = 1`; `zero_pow` twice gives `D_b(0) = 1`; `intermediate_value_Ioo`
on `Icc (−1) 0` puts `0 ∈ Ioo (−1) 1` in the image of `Ioo (−1) 0`.

`lake build` **exit 0** (3203 jobs). `lake test` **exit 0** (3201 jobs).

## What is **not** formalised

Stated plainly, because the Lean file must not be read as more than it is:

1. **`D_{a,b} = Y^{(a,b)}_{(a,b)}`.** On paper this is the `m = b` case of Theorem C. In Lean
   `t^b − t^{b−1} + 1` is simply *written down*. `Y`, Hall–Littlewood `P`, Kostka–Foulkes and
   charge **have no Lean definitions in this project**. A paper↔Lean dictionary is not a Lean
   definition.
2. **Root ⇒ non-cyclotomic.** Needs "every root of a cyclotomic polynomial has modulus 1". This
   was the session's stretch goal and I did not start it; the must-have plus the instrument work
   filled the hour. `exists_root_Ioo` is the input it would consume.
3. **The product-form null itself** — the quantification over all `λ` with `ℓ(λ) ≤ 2` and all
   two-part `ρ`.

`unproved ≠ unformalised`. The registry parent `thm-D-product-form-obstruction` therefore stays
`proved`, **not** `lean-verified`; only the new child is `lean-verified`.

Carried forward from the paper unchanged: Theorem D's conclusion needs **one** non-cyclotomic
witness and the odd-`b` family supplies infinitely many, so `gap:evenb` (even `b`, *computed*
for `3 ≤ b ≤ 15`, not proved) does not weaken the citation Rick asked for. Nothing in the Lean
file implies otherwise.

## Instruments — both arms, every reading

| instrument | clean arm | planted arm | alive? |
|---|---|---|---|
| `#print axioms` → `sorryAx`, 5 decls | **0** | **2** | **yes** |
| comment-stripped `sorry` grep | **0** | — (raw grep also 0; no prose "sorry" in this file) | n/a |
| `grep "declaration uses 'sorry'"` | 0 | **0 with a sorry planted** | **DEAD** |
| `lake build` exit code | 0 | **0 with a sorry planted** | **DEAD for sorries** |
| `lake test` exit code | 0 | **1** (one `#guard` value 3→4) | **yes** |
| `trustcheck … validate` | OK | **see below** | **partly** |

The planted `sorry` went in `D_eval_neg_one`; `sorryAx` appeared **there and in
`exists_root_Ioo`**, i.e. it propagates to the dependent — that propagation is the reason this
instrument is trustworthy and the build exit code is not. Restored reading: 0.

Confirming 10-07's finding rather than re-discovering it: the standard
`declaration uses 'sorry'` grep is **still a constant function in this toolchain**, and
`lake build` exits **0** with a `sorry` present. Neither can be used.

### Final `#print axioms`

```
'TworowD4Kernel.D'                        depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.continuous_D'             depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.D_eval_neg_one'           depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.D_eval_zero'              depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.exists_root_Ioo'          depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.D_four_eval_neg_one'      depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.D_four_pos_of_mem_Icc'    depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.not_exists_root_Ioo_four' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Standard three on all 8. `sorryAx` 0, `Lean.ofReduceBool`/`native_decide` 0.

## Ablation — and it ablated

**Mechanical arm.** Dropped `(hodd : Odd b)` from the statement and deleted its single use.
Pattern asserted present *before* mutating (both counts printed as `1`); **substitutions = 2**;
`lake build` **exit 1**:

```
error: TworowD4Kernel/NonCyclotomicRoot.lean:100:50: unsolved goals
b : ℕ
hb : 3 ≤ b
hcont : ContinuousOn (D b) (Icc (-1) 0)
⊢ Odd b
```

The residual goal is literally `⊢ Odd b` — the hypothesis, and nothing else, is what the proof
cannot do without.

**Mathematical arm, which is the stronger one.** A build failure only shows *my proof* needs the
hypothesis. `not_exists_root_Ioo_four` shows the **statement with `Odd b` deleted is false**:
`b = 4` satisfies `3 ≤ b`, and `t^4 − t^3 + 1` has **no** root in `(−1,0)` because
`D_four_pos_of_mem_Icc` proves it is strictly positive on all of `[−1,0]`. So the hypothesis is
not an artefact of the proof route — it is necessary for the theorem.

Executable shadows in `TworowD4KernelTests.lean` (over `ℤ`, where `decide` applies): the two
signs at `b = 3` (`−1 < 0`, `+1` at `t = 0`) and the reversed sign at `b = 4` (`+3 > 0`).

## A finding about trustcheck: it does not resolve the `lean` field

Two arms on the validator, because I was about to report a `lean-verified` grade on its say-so.

- **Arm A** — replaced `"lean": "TworowD4Kernel.exists_root_Ioo"` with
  `"TworowD4Kernel.no_such_declaration_xyz"`. Reading: **`OK: … is valid`.**
- **Arm B** — misspelled the trust value as `lean-verfied`. Reading: **`2 problem(s)`** — the
  invalid enum *and* the boundary rule on the parent.

So the validator is alive on the trust enum and on parent/child boundary rules, and **blind on
`lean`**: it never checks that the named declaration exists, let alone that it is the statement
claimed. *The green on a `lean-verified` node rests entirely on my own `#print axioms` reading,
not on the validator.* `trustcheck` grades the shape of the tree, not the Lean. Worth knowing
before the next session reads an `OK` as corroboration of a `lean` field.

Note also, confirming the brief: `--files-dir .`, not `proofs`, which double-prefixes.

## Honest friction, for the record

- The pre-flight "assert the pattern is present before mutating" guard **fired on me**: I
  predicted `rw [D, h1, h2]` occurred once; it occurs **twice** (`D_eval_neg_one` *and*
  `D_eval_zero`), and the mutation would have hit the wrong one or both. Re-anchored on a
  unique string. The guard earned its place in the protocol on its first use this session.
- `lake` is **not on `PATH`** and `~/.elan` does not exist; the binary is at
  `/home/clio/projects/.elan/bin/lake` (`ELAN_HOME=/home/clio/projects/.elan`). My first build
  command read `LAKE_EXIT=1` — which was `which lake` failing, not Lean. Exactly the
  exit-code-grades-the-command-that-ran trap, caught because the brief told me to expect it.
- `import Mathlib.Analysis.SpecialFunctions.Polynomials` **does not exist** in Mathlib v4.30.0
  (`bad import`); the needed imports are `Mathlib.Topology.Order.IntermediateValue`,
  `Mathlib.Topology.Instances.Real.Lemmas`, `Mathlib.Topology.Algebra.Ring.Basic`.

## Registry

Added `two-part-green-polynomials.json` → `thm-D-product-form-obstruction` →
**`lean-root-in-Ioo-odd-b`**, `trust: lean-verified`,
`lean: TworowD4Kernel.exists_root_Ioo`,
`file: lean/tworow_d4_kernel/TworowD4Kernel/NonCyclotomicRoot.lean`. Parent unchanged at
`proved`. Validator `OK` — with the caveat above about what that `OK` does and does not test.
Backup kept at `…json.bak-1007c2-lean`.
