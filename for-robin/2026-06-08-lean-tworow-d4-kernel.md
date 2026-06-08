# Lean session 2026-06-08 — arithmetic kernel of the b ≡ 1 (mod 4) two-row d=4 proof

Robin —

First Lean session with a working toolchain. There was no Lean/elan/lake installed in the
container and no prior Mathlib cache, so I bootstrapped from scratch: installed `elan`, set up a
fresh `lake` project against **Mathlib v4.30.0**, pulled the precompiled oleans
(`lake exe cache get`), and formalised the two arithmetic lemmas LEAN.md asked for. Both compile
with **zero errors, zero warnings, no `sorry`**, and `#print axioms` shows only the three standard
Mathlib axioms (`propext`, `Classical.choice`, `Quot.sound`).

## What's proved

**PRIMARY (§2.1, "the j=1 term is exact").** The falling-factorial / binomial bridge:
```lean
theorem descFactorial_eq_factorial_mul_self_mul_choose_pred (m R : ℕ) (h : R + 1 ≤ m) :
    Nat.descFactorial m (R + 1) = R ! * (m * Nat.choose (m - 1) R)
```
i.e. `(m)_{R+1} = R! · (m · C(m−1, R))`. This is what makes `τ₁(m) = (m)_{R+1}/R!` an integer
and pins `v₂(τ₁) = v₂((m)_{R+1}) − v₂(R!)`. Proof: write `m = k+1`, apply
`Nat.descFactorial_eq_factorial_mul_choose`, `Nat.factorial_succ`, `Nat.add_one_mul_choose_eq`,
close with `ring`. Clean.

**SECONDARY (§2.3, "the 2-adic doubling identity").**
```lean
theorem padicValNat_two_factorial_two_mul (h : ℕ) :
    padicValNat 2 (2 * h)! = h + padicValNat 2 h !
```
i.e. `v₂((2h)!) = h + v₂(h!)`. Pleasant surprise: Mathlib already has the *multiplicative* form
of Legendre's theorem, `Nat.padicValNat_factorial_mul : padicValNat p (p*n)! = padicValNat p n! + n`,
so I did **not** need the digit-sum route LEAN.md anticipated (`sub_one_mul_padicValNat_factorial`
+ `Nat.digits`). The whole proof is one rewrite plus `add_comm`. The digit-sum lemmas are there
too (`sub_one_mul_padicValNat_factorial`) if a future, less-aligned valuation identity needs them.

## Where it lives

- Compiling project: `~/projects/lean/tworow_d4_kernel/` (full lake project, has the `.lake`
  build tree — local only).
- Standalone source copy: `~/projects/lean/2026-06-08-tworow-d4-kernel.lean`.
- Pushed to **clio-vega/proofs**: `2026-06-08-lean-tworow-d4-kernel/` (source + lakefile +
  lean-toolchain + lake-manifest + README with build instructions; `.lake` omitted). If you open
  it in VS Code, `lake exe cache get` then `lake build` reproduces it in ~30 s.

## Honest status

Nothing is left as `sorry`. The scope is deliberately just the two arithmetic lemmas — **not** the
`I_b(m)` generating function or the term-domination argument (§2.4), which need machinery Mathlib
doesn't carry and would be a multi-session effort. These two are the load-bearing valuation facts,
and they are now machine-checked. If you'd like, the natural next Lean target is the closed form
for `D_j` in §2.3 (which *uses* the doubling identity) — that's a discrete bookkeeping lemma that
might be in reach once I have the surrounding definitions stated.

— Clio
