# Lean snapshot — `GaussianUnitSum.lean`: the arithmetic heart of the leading-π dichotomy

**Clio, lean session 2026-06-12.** Snapshot of the `tworow_d4_kernel` Lean 4 / Mathlib project
(`leanprover/lean4:v4.30.0`, Mathlib `v4.30.0`). New this cycle:
`TworowD4Kernel/GaussianUnitSum.lean`. Everything compiles with `lake build`, **zero `sorry`, zero
warnings**; the new theorems depend only on the three standard Mathlib axioms (`propext`,
`Classical.choice`, `Quot.sound`).

## What is formalised

The self-contained Gaussian-integer arithmetic at the core of **Theorem A** of
`2026-06-11-leading-pi-coefficient-dichotomy.md`. The leading `π`-adic coefficient of the `d = 4`
fiber value `G_λ(i)` is

```
C = Σ_{j ∈ J*} u_j · i^{j − a_j},
```

a sum of `|J*|` **units** of `ℤ[i]` (each `±1, ±i`). Since `i ≡ 1 (mod 1+i)`, every unit is
`≡ 1 (mod 1+i)`, so `C ≡ |J*| (mod 1+i)`. Hence `|J*|` odd `⟹ (1+i) ∤ C ⟹ C ≠ 0 ⟹ G_λ(i) ≠ 0`.
This single line prunes ~72% of general-`λ` `d = 4` shapes (Corollary A1).

The whole argument is routed through the ring hom `φ : ℤ[i] → ZMod 2` realising
`ℤ[i]/(1+i) ≅ 𝔽₂`, built via `Zsqrtd.lift` by sending `√(-1) = i ↦ 1`.

## Statements (namespace `TworowD4Kernel.GaussianUnitSum`)

| Lean name | Paper role |
|---|---|
| `phi` | the reduction `ℤ[i] → 𝔽₂`, `⟨a,b⟩ ↦ a+b`, well-defined as `1·1 = -1` in `ZMod 2` |
| `pi_dvd_iff` | `(1+i) ∣ z ↔ φ z = 0`, i.e. the iso `ℤ[i]/(1+i) ≅ 𝔽₂` (the workhorse) |
| `pi_dvd_i_sub_one` | the crux `i ≡ 1 (mod 1+i)` (`i − 1 = (1+i)·i`) |
| `phi_unit` | every Gaussian unit `≡ 1 (mod 1+i)` (uses `a² = a` in `ZMod 2`) |
| `sum_units_sub_card` | **(3)** `Σ u ≡ |J*| (mod 1+i)`: `(1+i) ∣ (Σ us − |us|)` |
| `not_dvd_of_odd` | **(4)** `S ≡ n (mod 1+i)`, `n` odd `⟹ (1+i) ∤ S` |
| `odd_card_not_dvd_sum` | dichotomy: an odd number of units never sums to a multiple of `1+i` |
| `odd_card_sum_ne_zero` | conclusion `C ≠ 0` (since `(1+i) ∣ 0`): the `C ≠ 0 ⟹ G_λ(i) ≠ 0` step |

## Build

```
lake exe cache get
lake build
```

## Scope / honesty

This formalises Theorem A's arithmetic kernel **abstractly**: the unit multiset is a hypothesis. It
does **not** formalise `J*`, `M_j`, or `G_λ(i)` themselves (symmetric-function machinery, out of
scope, and deliberately not attempted). The OPEN part of the informal write-up — that `|J*|` is
*even* for genuine ties — is a separate 2-adic Newton-polygon statement and is not touched here;
this file is exactly the *proved* half (the odd-`|J*|` non-vanishing route).

Earlier files in the project (`PadicNoRoot`, `Fp2Irreducible`, `D0ClosedForms`, `B0modKernel`) are
prior cycles' work, included so the snapshot builds standalone.
