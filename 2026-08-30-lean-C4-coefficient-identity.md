# LEAN 2026-08-30 — the algebraic core of C4, machine-checked

**Project:** `~/projects/lean/tworow_d4_kernel/` (Lean 4 `v4.30.0`, Mathlib `v4.30.0`)
**New file:** `TworowD4Kernel/QuantumInteger.lean` (164 lines), namespace `Clio.QuantumInteger`
**Paper proof:** `~/projects/proofs/2026-08-14-C4-iijima-B1.tex`, §3 (Theorem 3.1 and its proof)
**Registry:** `proofs/registry/fock-ribbon-sign-operator.json`, new node
`root/C4-explicit-quotient/C4-coefficient-identity`, `trust: lean-verified`

## Target declaration

```lean
theorem Clio.QuantumInteger.C4_coefficient_identity
    {R : Type*} [CommRing R] (q : Rˣ) (n : ℕ) :
    (-qpow q 1) ^ n - (-qpow q (-1)) ^ n
      = (-1) ^ n * ((qpow q 1 - qpow q (-1)) * qInt q n)
```

with

```lean
def qpow (q : Rˣ) (n : ℤ) : R := ((q ^ n : Rˣ) : R)
def qInt (q : Rˣ) (n : ℕ) : R := ∑ k ∈ Finset.range n, qpow q ((n : ℤ) - 1 - 2 * k)
```

`q` is a **unit** of an arbitrary commutative ring, so `q⁻¹` exists without inverting anything
else in `R` and negative powers are `zpow`. The intended instance is `R = ℤ[q, q⁻¹]`, the
coefficient ring of level-1 Uglov Fock space `𝓕_e`.

## What builds, sorry-free

Everything. `lake build TworowD4Kernel.QuantumInteger` succeeds with **zero sorries and zero
warnings**. Eight declarations:

| declaration | statement |
|---|---|
| `qpow_zero` | `q^0 = 1` |
| `qpow_add` | `q^a · q^b = q^{a+b}` |
| `qpow_natPow` | `(q^a)^n = q^{an}` (`n : ℕ`) |
| `qpow_inv` | `(q⁻¹)^a = q^{-a}` |
| `qInt_zero`, `qInt_one_index` | `[0]_q = 0`, `[1]_q = 1` |
| **`qpow_sub_qpow_inv`** | **`q^n - q^{-n} = (q - q⁻¹) · [n]_q`** |
| **`C4_coefficient_identity`** | **`(-q)^n - (-q⁻¹)^n = (-1)^n ((q - q⁻¹) [n]_q)`** |
| `qInt_one` | `[n]_q |_{q=1} = n` |
| `qInt_inv` | `[n]_{q⁻¹} = [n]_q` (bar-invariance) |

`qpow_sub_qpow_inv` is the only real content, exactly as anticipated. The proof is the telescope
`(q - q⁻¹) · q^{n-1-2k} = q^{n-2k} - q^{n-2(k+1)}` summed over `k ∈ range n`, closed by
`Finset.sum_range_sub'`. The friction was where LEAN.md predicted — the `ℤ`-exponent arithmetic
`n - 1 - 2k` and its `ℕ → ℤ` casts, not the algebra. Two build errors total, both cast/rewrite
bookkeeping (`neg_pow` matching `(-1)^n` as well as `(-q)^n`; a `congr 1` closing a goal early).

`qInt_one` and `qInt_inv` are the two sanity checks LEAN.md asked for as the honest bonus. They
pin the exponent convention `n-1-2k` against a misreading — the specific failure mode this
programme keeps hitting.

## `#print axioms`

```
'Clio.QuantumInteger.qpow_sub_qpow_inv'     depends on axioms: [propext, Classical.choice, Quot.sound]
'Clio.QuantumInteger.C4_coefficient_identity' depends on axioms: [propext, Classical.choice, Quot.sound]
'Clio.QuantumInteger.qInt_one'              depends on axioms: [propext, Classical.choice, Quot.sound]
'Clio.QuantumInteger.qInt_inv'              depends on axioms: [propext, Classical.choice, Quot.sound]
'Clio.QuantumInteger.qInt_one_index'        depends on axioms: [propext, Classical.choice, Quot.sound]
```

The standard three, nothing else.

## What this does *not* verify — the honest boundary

C4 as stated in the paper is an **operator identity on `𝓕_e`**:

> `P_e = B_{-1}^{[e]} + (q - q⁻¹) C_e^{(1)}`.

Lean holds the coefficient identity only. Three steps of the paper proof stay informal, and they
are exactly the steps that were never in scope here:

1. that the `|μ⟩`-coefficient of `P_e|λ⟩` is `(-q)^{h(μ/λ)}` and of `B_{-1}^{[e]}|λ⟩` is
   `(-q⁻¹)^{h(μ/λ)}` (Definitions 2.x, and Iijima arXiv:1207.6161 eq. (9) at `m = -1`);
2. that `μ` ranges over `λ + e`-ribbon and `h` is the ribbon height — no ribbons, partitions or
   abaci exist in the Lean file;
3. the `ℤ[q,q⁻¹]`-linear extension from coefficients to operators.

So the registry node is a **child** of `C4-explicit-quotient`, not a re-grading of it: the parent
stays `peer-reviewed` (Lyra, 2026-08-13). Marking the parent `lean-verified` would be trust
inflation of precisely the kind the registry README warns about.

## Finding: the project has an import cycle, and C1–C4 boundary cannot be compiled

LEAN.md asked me to spend five minutes answering Rick's twice-made request — are the `c = 1..4`
2-adic content lemmas axiom-free? — with `#print axioms` on the four boundary theorems.

**I cannot answer it, and the reason is a real defect, not a missing five minutes.**

```
TworowD4Kernel                        (root aggregator)
  → TworowD4Kernel.ThreeRowC1Boundary
  → TworowD4Kernel.HookKummerLemmas
  → TworowD4Kernel.D0ClosedForms
  → TworowD4Kernel                    ← cycle
```

`D0ClosedForms.lean:6` is `import TworowD4Kernel`, and the root file (last touched 2026-07-04) now
imports the four boundary modules. Under the shipped `globs = ["TworowD4Kernel.+"]` the root is
not a library target at all, so `D0ClosedForms` fails with *"object file `TworowD4Kernel.olean`
does not exist"* and takes `HookKummerLemmas`, `NumberLemmaC2`, `CompensationLemma`,
`SubsetIdentityGeneralC` and all four `ThreeRowC*Boundary` modules down with it. Adding the root
to the globs (tried, then reverted) does not help — it turns the missing-target error into an
honest `build cycle detected`.

Consequences, stated plainly:

* **`lake build` on this project currently fails.** It was described in LEAN.md as sorry-free and
  building; the sorry-free part I confirmed by inspection (`grep` finds no `sorry`, no `axiom`, no
  `native_decide` anywhere in the sources), but *nothing above `D0ClosedForms` compiles today*.
* The four boundary theorems' axiom dependencies are therefore **unknown**, and their `sorry`-free
  advertisement in their own module docstrings is currently **unchecked by the compiler**.
* The repair is not a one-liner: `D0ClosedForms` genuinely needs the two lemmas that live in the
  root file (`descFactorial_eq_factorial_mul_self_mul_choose_pred`,
  `padicValNat_two_factorial_two_mul`). Breaking the cycle means moving those two lemmas out of
  the root aggregator into a base module — a refactor, and out of scope for a session whose rule
  was "stay with the target".

This outranks the axioms answer in importance, so it is recorded here rather than patched.

## Session infrastructure note

**The container had no Lean toolchain at all.** No `elan`, no `lake`, no `lean`, no Mathlib
checkout, no `.olean` anywhere on the filesystem — presumably lost with the 18-day outage. About
half the session went to reinstalling: `elan` + Lean `v4.30.0` (aarch64), Mathlib at the pinned
rev, and `lake exe cache get` (8459 files). One trap worth recording: **`/home/clio/.cache` is
owned by `root` and not writable by `clio`**, so `cache get` dies with `permission denied (error
code: 13)`. Workaround, needed every session until the container image is fixed:

```
export XDG_CACHE_HOME=/home/clio/.lean-cache
```

Cold-cache elaboration of one small file against `import Mathlib.Tactic` is ~250–570 s, so budget
two or three build attempts per hour, not fifty.
