# Lean session 2026-06-09 — arithmetic kernel of the b ≡ 0 (mod 4) two-row d=4 proof

Robin —

This is the file LEAN.md asked for last cycle and that didn't get delivered: `B0modKernel.lean`,
the arithmetic kernel of the `b ≡ 0 (mod 4)` proof
(`2026-06-08-tworow-d4-b0mod4-proved.md`). It's now in the existing `tworow_d4_kernel` project,
imported by the root module, building alongside the `b ≡ 1` lemmas. **Zero errors, zero warnings,
no `sorry`**; `#print axioms` shows only the three standard axioms (`propext`, `Classical.choice`,
`Quot.sound`) on all three new declarations.

## What's proved

**PRIMARY (§3, "the `h` consecutive integers ÷ `h!` engine").**
```lean
theorem padicValNat_two_choose_eq (R h : ℕ) (hh : h ≤ R) :
    padicValNat 2 (Nat.choose R h) =
      padicValNat 2 (Nat.descFactorial R h) - padicValNat 2 (Nat.factorial h)
```
i.e. `v₂(C(R, h)) = v₂((R)_h) − v₂(h!)`. This is the valuation identity sitting under every closed
form `D₀(j)` in §3 (Lemma 1) — a product of `h` consecutive integers is divisible by `h!`, and the
quotient is `C(R, h)`. Proof is three rewrites: `Nat.descFactorial_eq_factorial_mul_choose` turns
`(R)_h` into `h! · C(R, h)`, then `padicValNat.mul` (needs `Fact (Nat.Prime 2)` and both factors
nonzero — `h ≤ R` is exactly what makes `C(R, h) ≠ 0`), then `Nat.add_sub_cancel_left`. Clean and
short.

**SECONDARY (§6, "three odds sum to an odd").** The parity count that rescues the law when strict
domination fails (for `m` odd the minimal-valuation terms are `{1, 2, 3}`, three odds):
```lean
theorem odd_add_three {a b c : ℕ} (ha : Odd a) (hb : Odd b) (hc : Odd c) : Odd (a + b + c)
theorem odd_sum_of_odd_card {ι : Type*} (s : Finset ι) (f : ι → ℕ)
    (hf : ∀ i ∈ s, Odd (f i)) (hcard : Odd s.card) : Odd (∑ i ∈ s, f i)
```
The first is the bare three-term instance actually invoked. The second is the general packaging
("an odd-cardinality sum of odds is odd"), proved by reducing mod 2: each `f i ≡ 1`, so
`∑ f i ≡ s.card ≡ 1 (mod 2)` via `Finset.sum_nat_mod` + `Finset.sum_const`.

## Where it lives

- Compiling project: `~/projects/lean/tworow_d4_kernel/` (full lake project, local only).
- Standalone source copy: `~/projects/lean/2026-06-09-tworow-d4-b0mod-kernel.lean`.
- Pushed to **clio-vega/proofs**: `2026-06-09-lean-tworow-d4-b0mod-kernel/` (full updated project
  snapshot: root module + the new `TworowD4Kernel/B0modKernel.lean` + lakefile/toolchain/manifest +
  README; `.lake` omitted). `lake exe cache get` then `lake build` reproduces it.

## Honest status / scope

Nothing is `sorry`. As before, the scope is deliberately the arithmetic kernel only — **not** the
`I_b(m)` generating function, the term decomposition `v₂(τ_j) − a = D₀(j) + S_j(m)`, or the tie
classification (Lemma 3). Those need the surrounding combinatorial definitions Mathlib doesn't
carry. With the `b ≡ 1` doubling identity (last cycle) and this `b ≡ 0` choose/descFactorial
engine, the two load-bearing valuation facts of *both* proved residue classes are now
machine-checked.

A natural next Lean step, if you want to push deeper: state `D₀(j)` for `j` odd/even (Lemma 1) and
prove the two closed forms `D₀ = v₂(C(R, h))` / `D₀ = v₂(C(R, h−1)) − v₂(h)` — both are now within
reach since they are assembled from exactly the doubling identity and this choose-valuation lemma
plus the `β_j = R − h` / `R + 1 − h` index bookkeeping.

— Clio
