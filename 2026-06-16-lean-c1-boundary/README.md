# Lean snapshot — three-row `c = 1` boundary lemma, end-to-end

**Date:** 2026-06-16/17 (lean session) · **Project:** `lean/tworow_d4_kernel`, Mathlib v4.30.0

`ThreeRowC1Boundary.lean` assembles the boundary lemma for the three-row tie family `λ = (a, b, 1)`
into a single `sorry`-free declaration `threerow_c1_boundary`, closing Gap 3 (the boundary residual)
of the `c = 1` interior theorem.

## Status
- `lake build`: 0 errors, 0 warnings.
- `#print axioms threerow_c1_boundary` → `[propext, Classical.choice, Quot.sound]` only (no `sorry`,
  no extra axioms). Same for the two parity sub-theorems and both digit-sum helpers.

## What is machine-checked here (pure 2-adic number theory)
- `two_mul_sum_digits_le` — **Lemma A**, the digit-sum envelope `2 s₂(k) ≤ k − 1` (`k ≥ 4`), strong
  induction on the binary representation. This is the engine that closes the `c = 1` boundary — the
  brief's expectation that Lemma F drives `c = 1` was off; §2 of the boundary note actually uses the
  digit-sum envelope. (Lemma F drives `c = 2, 3`.)
- `sum_digits_two_add_le` — digit-sum subadditivity `s₂(x+y) ≤ s₂(x) + s₂(y)` (= the carry count
  `γ ≥ 0`), strong induction; the both-odd case folds in the `+1` instance.
- `threerow_c1_boundary_aeven` / `_aodd` — `Δ ≥ 2` / `Δ ≥ −1`, the parity closures (the `b = 3`
  corner is handled by hand: `m` odd ⟹ `4 ∣ a+2`; `m` even ⟹ the predecessor bound `s₂(k) ≤
  s₂(k−1)+1`).
- `threerow_c1_boundary` — the assembled statement `V < val(b+1)` under the threshold `V = val(0)−θ`.

## What remains MN-verified (NOT in Lean, by design)
The Proposition-1 closed form `Δ = (b−3) + 2 v₂(a+2) + 2(s₂(m−b) − s₂(a))` is taken as the hypothesis
`hΔ`: it rests on the Murnaghan–Nakayama closed forms for `M_0 = f^λ` (hook formula) and
`M_{b+1} = b`, symmetric-function facts out of Mathlib's reach — same scoping as every prior session.
