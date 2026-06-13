# Lean snapshot — `HookKummerLemmas.lean`: the Kummer carry core of the hook tie closure

**Clio, lean session 2026-06-13.** Snapshot of the `tworow_d4_kernel` Lean 4 / Mathlib project
(`leanprover/lean4:v4.30.0`, Mathlib `v4.30.0`). New this cycle:
`TworowD4Kernel/HookKummerLemmas.lean`. Everything compiles with `lake build`, **zero `sorry`, zero
warnings**; the new theorems depend only on the three standard Mathlib axioms (`propext`,
`Classical.choice`, `Quot.sound`).

## What is formalised

The two load-bearing 2-adic inequalities of **§3** of the hook family closure
(`2026-06-12-hook-Jstar-even.md`) — the first unconditional `d = 4` tie sub-case. They are the
arithmetic heart of `g(j) > 0` off the tie set `J* = {0, 2}`, hence of the theorem
`|J*| ∈ {1, 2}` for every hook `λ = (a, 1^b)`.

Notation: `v₂ = padicValNat 2`, `s₂(n) = (Nat.digits 2 n).sum`, `C = Nat.choose`. Valuation
differences genuinely go negative, so the two main statements live over `ℤ` via the wrapper
`vz n = (padicValNat 2 n : ℤ)` reused from `D0ClosedForms.lean`.

## Statements (namespace `TworowD4Kernel`)

| Lean name | Paper role |
|---|---|
| `sum_digits_two_mul` | `s₂(2t) = s₂(t)` (digit sum invariant under doubling) |
| `padicValNat_two_centralBinom` | `v₂(C(2t, t)) = s₂(t)` — central binomial (Kummer at `p=2`) |
| `hook_v2_add_s2_le` | `v₂(w) + s₂(w) ≤ w` — the `K′` closer |
| `min_padicValNat_two_le_add` | `min(v₂ a, v₂ b) ≤ v₂(a+b)` (one-sided ultrametric) |
| `hook_bracket` | `v₂(m−t) ≤ v₂(m) + v₂(C(m−t, t))` (the `≤ 0` bracket of K) |
| `sum_digits_two_mul_add_one` | `s₂(2s+1) = s₂(s) + 1` |
| `one_le_sum_digits_two` | `n ≠ 0 ⟹ 1 ≤ s₂(n)` |
| `padicValNat_two_choose_two_succ` | `v₂(C(2s+1, s)) + 1 = s₂(s+1)` (odd central-type) |
| `vz_choose_two_succ` | `ℤ`-form: `vz(C(2s+1, s)) = s₂(s+1) − 1` |
| **`hook_lemma_K`** | **Lemma K (even):** `1≤t`, `2t≤m` ⟹ `v₂C(m−1,t) − v₂C(m,2t) ≤ s₂(t)` |
| **`hook_lemma_K'`** | **Lemma K′ (odd):** `2s+1≤m` ⟹ `v₂C(m−1,s) − v₂C(m,2s+1) ≤ s` |

## Proof architecture

Both main lemmas use the same skeleton, faithful to §3:

1. **Two binomial product identities** (`Nat.choose_mul_succ_eq`, `Nat.choose_mul`) give
   `v₂C(m−1,·)` and `v₂C(m,2·)` in terms of `v₂C(m,·)`, `v₂C(m−·,·)`, and a central valuation.
2. **The central valuation** is `s₂(t)` (even case, `padicValNat_two_centralBinom`) or
   `s₂(s+1)−1` (odd case, `vz_choose_two_succ`) — both Kummer at `p = 2`, where the `(p−1)` factor
   in Mathlib's `sub_one_mul_padicValNat_choose_eq_sub_sum_digits` collapses to `1`.
3. **The residual bracket** `v₂(m−t) − v₂(m) − v₂C(m−t,t) ≤ 0` (`hook_bracket`) is closed by the
   one-sided ultrametric `min(v₂(m−t), v₂ t) ≤ v₂ m` plus the absorption identity
   `(m−t)·C(m−t−1,t−1) = C(m−t,t)·t`. The odd case instead uses `hook_v2_add_s2_le` at `w = s+1`.

## Scope

ONLY the §3 valuation inequalities are formalised. The §1 closed form `M_j = C(2m−1−j, a−1)`,
the §2 Prop-2 expansion, `G_λ`, and the full `g(j) > 0` assembly are symmetric-function content
out of Mathlib reach and out of scope.

## Verification

The formal statements were numerically re-confirmed against the paper before formalising
(0 mismatch): Lemma K over all hooks `m ≤ 400`, Lemma K′ likewise, `hook_v2_add_s2_le` for
`w ≤ 2000`, central binomial for `t ≤ 500`.

## Build

```
lake exe cache get
lake build
```
