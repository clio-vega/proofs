# Lean 4 formalisation — arithmetic kernel of the two-row d=4 law (b ≡ 1 mod 4)

**Date:** 2026-06-08 (lean session)
**Source proof:** `../2026-06-07-tworow-d4-b1mod4-proved.md` (§2.1 and §2.3)
**Toolchain:** Lean 4 `v4.30.0`, Mathlib `v4.30.0` (see `lean-toolchain`, `lake-manifest.json`)

## What is formalised

Two self-contained arithmetic lemmas that are load-bearing in the 2-adic proof of the two-row
d=4 fiber-vanishing law `G_{(2m−b,b)}(i) = 0 ⟺ (2,2)` for the infinite family `b ≡ 1 (mod 4)`.
The full theorem rests on generating-function machinery that is not in Mathlib; these two lemmas
are exactly the valuation bookkeeping a skeptic would most want machine-checked. **Both are fully
proved — no `sorry`, no `axiom` beyond Mathlib's standard `{propext, Classical.choice, Quot.sound}`.**

### 1. `descFactorial_eq_factorial_mul_self_mul_choose_pred` (§2.1, "the j=1 term is exact")

```lean
theorem descFactorial_eq_factorial_mul_self_mul_choose_pred (m R : ℕ) (h : R + 1 ≤ m) :
    Nat.descFactorial m (R + 1) = R ! * (m * Nat.choose (m - 1) R)
```

The falling-factorial / binomial bridge `(m)_{R+1} = R! · (m · C(m−1, R))`. This is what makes
`τ₁(m) = (m)_{R+1}/R!` an honest integer and pins `v₂(τ₁(m)) = v₂((m)_{R+1}) − v₂(R!)`.

*Proof:* write `m = k + 1` (forced by `R + 1 ≤ m`), then `descFactorial m (R+1) = (R+1)!·C(m,R+1)`
(`Nat.descFactorial_eq_factorial_mul_choose`), `(R+1)! = (R+1)·R!` (`Nat.factorial_succ`), and
`(k+1)·C(k,R) = C(k+1,R+1)·(R+1)` (`Nat.add_one_mul_choose_eq`); `ring` finishes.

### 2. `padicValNat_two_factorial_two_mul` (§2.3, "the 2-adic doubling identity")

```lean
theorem padicValNat_two_factorial_two_mul (h : ℕ) :
    padicValNat 2 (2 * h)! = h + padicValNat 2 h !
```

The valuation step `v₂((2h)!) = h + v₂(h!)`. It is the `p = 2` instance of the Legendre-theorem
corollary `Nat.padicValNat_factorial_mul : padicValNat p (p*n)! = padicValNat p n! + n`, so the
proof is a one-line rewrite plus `add_comm`. (No digit-sum manipulation was needed — Mathlib
already carries the multiplicative form of Legendre's theorem.)

## How to build

This directory is a slim, reproducible snapshot (source + config only; the multi-GB `.lake`
build tree is omitted). To check it:

```sh
elan toolchain install $(cat lean-toolchain)        # if not present
lake exe cache get                                   # fetch precompiled Mathlib oleans
lake build                                           # ~30 s after cache
```

You will need a top-level `TworowD4Kernel.lean` import target; the lakefile's default target
already points at it. Place `TworowD4Kernel.lean` at the project root (as here).
