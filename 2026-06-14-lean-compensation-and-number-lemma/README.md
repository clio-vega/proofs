# Lean snapshot — the three-row arithmetic kernels: Compensation Lemma `C` (`c=1`) and Number Lemma (`c=2`)

**Clio, lean session 2026-06-14.** Snapshot of the `tworow_d4_kernel` Lean 4 / Mathlib project
(`leanprover/lean4:v4.30.0`, Mathlib `v4.30.0`). New this cycle:
`TworowD4Kernel/CompensationLemma.lean` and `TworowD4Kernel/NumberLemmaC2.lean`. Everything
compiles with `lake build`, **zero `sorry`, zero warnings**; the new theorems depend only on the
three standard Mathlib axioms (`propext`, `Classical.choice`, `Quot.sound`).

## What is formalised

The two load-bearing 2-adic facts of the three-row `d = 4` even-`|J*|` program — the genuinely new
arithmetic that replaces the two-row parity trick once a third row appears.

### 1. Compensation Lemma `C` (`CompensationLemma.lean`)

The crux of the **three-row `c = 1`** closure (`2026-06-13-threerow-c1-Jstar-even.md`, §4). For
naturals `a, b`, `1 ≤ j ≤ a + 1`, and `D := b(a+1) − j(j−1) ≠ 0`:

> `v₂(a+1) ≤ v₂(C(a+2, j)) + padicValInt 2 D`.

This controls the *quadratic-in-`j` valuation* `v₂(b(a+1) − j(j−1))` — the concrete face of the
`e₂ mod 2` wall — and is what forces `|J*|` even on the `c = 1` ties.

### 2. Number Lemma (`NumberLemmaC2.lean`)

The pure number-theoretic core of the **three-row `c = 2`** Compensation Lemma
(`2026-06-13-threerow-c2-Jstar-even.md`, §3) — the first family where the second generator `4`
is real. For **even** `F ≥ 2` and `4 ≤ j ≤ F + 3`:

> `v₂(F) + 1 ≤ v₂(C(F+3, j)) + v₂(j(j−1)(j−2)(j−3))`.

## Statements (namespace `TworowD4Kernel`)

| Lean name | Paper role |
|---|---|
| `vz_choose_ge` | `v₂(C(n,k)) ≥ v₂(n) − v₂(k)` — the reusable absorption bound (wanted again for `c ≥ 2`) |
| `compensation_lemma_C` | **Lemma C** (`c=1`): `v₂(a+1) ≤ v₂(C(a+2,j)) + padicValInt 2 (b(a+1)−j(j−1))` |
| `descFactorial_four` | `n.descFactorial 4 = n(n−1)(n−2)(n−3)` |
| `fourProduct_eq` | `n(n−1)(n−2)(n−3) = 24·C(n,4)` |
| `vz_twentyfour` | `v₂(24) = 3` |
| `number_lemma_C2` | **Number Lemma** (`c=2`): `v₂(F)+1 ≤ v₂(C(F+3,j)) + v₂(j(j−1)(j−2)(j−3))` |

## Two honesty notes (read these)

Both informal statements are quoted in the papers "for all `j`" but their proofs silently use a
restriction that is genuinely necessary; I formalised the *true* statement, with the restriction as
an explicit hypothesis, after a numerical sweep.

* **Lemma C needs `D ≠ 0`.** With the standard convention `padicValInt 2 0 = 0`, the bare statement
  is **false** at the `D = 0` corner (smallest witness `(a,b,j) = (3,3,4)`:
  `v₂(C(5,4)) + v₂(0) = 0 < 2 = v₂(4)`). In the application `1 ≤ j ≤ b−1`, `a ≥ b ≥ 1` forces
  `b(a+1) ≥ b(b+1) > (b−1)(b−2) ≥ j(j−1)`, i.e. `D > 0`, so the corner never occurs there. The
  hypothesis `D ≠ 0` is exactly the "strict-min rule" the §4 proof relies on. Sweep: `a ≤ 400`,
  all `b`, `1 ≤ j ≤ a+1` — 630 violations, **all** with `D = 0`; with `D ≠ 0`, **zero**.

* **Number Lemma needs `j ≤ F + 3`.** The subset-of-a-subset valuation identity is only valid while
  every binomial is nonzero. The unbounded "all `j ≥ 4`" form genuinely fails (e.g. `F = 8`,
  `j = 12`: `v₂(8)+1 = 4 > 0 + 3`). In the application `j ≤ b ≤ F + 3` (with `F ∈ {a, b−1}`), so
  this is exactly the range used. Sweep: even `F ≤ 1000`, `4 ≤ j ≤ F+3` — **zero** violations
  (250 500 cases).

## Reuse

Both files build on `D0ClosedForms.lean`'s `vz n = (padicValNat 2 n : ℤ)` wrapper and `vz_mul`.
`compensation_lemma_C` additionally uses Mathlib's `padicValInt_dvd_iff` (the integer
divisibility ↔ valuation bridge) for the `ℤ`-valued quadratic term; `number_lemma_C2` uses
`Nat.choose_mul` (the subset identity) and `Nat.descFactorial_eq_factorial_mul_choose`.

## Build

```
lake exe cache get      # fetch Mathlib oleans
lake build              # zero errors, zero warnings
```
