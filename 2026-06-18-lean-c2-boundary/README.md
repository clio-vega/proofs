# Lean snapshot — three-row `c = 2` boundary lemma, end-to-end

**Date:** 2026-06-18 (lean session) · **Project:** `lean/tworow_d4_kernel`, Mathlib v4.30.0

`ThreeRowC2Boundary.lean` assembles the boundary lemma for the three-row tie family `λ = (a, b, 2)`
into a single `sorry`-free declaration `threerow_c2_boundary`, closing the boundary residual (Gap 3)
of the `c = 2` interior theorem. With `c = 1` (snapshot `2026-06-16-lean-c1-boundary/`) this makes
**both** complete three-row families end-to-end machine-checked.

## Status
- `lake build`: 0 errors, 0 warnings (full project: 2090 jobs, success).
- `#print axioms threerow_c2_boundary` → `[propext, Classical.choice, Quot.sound]` only (no `sorry`,
  no extra axioms). Same for the two per-index sub-theorems and the four helpers.

## What is machine-checked here (pure 2-adic number theory)
- `vz_factorial` — Legendre at `p = 2` in `ℤ`-form: `v₂(n!) = n − s₂(n)` (from Mathlib's
  `sub_one_mul_padicValNat_factorial`).
- `factorial_mul_prod_Icc` — `P! · ∏_{t=1}^{Q}(P+t) = (P+Q)!` (induction, `prod_Icc_succ_top`).
- `vz_prod_Icc` — **the bridge** `v₂(∏_{t=1}^{Q}(P+t)) = Q + s₂(P) − s₂(P+Q)`. This is the form in
  which the carry-sum + `v₂(Q!)` of §3 enters every `a`-even case; it lets `lemma_F` apply directly.
- `vz_dvd_le` — valuation monotonicity `x ∣ y ⟹ v₂(x) ≤ v₂(y)`. Closes the subtop `a`-even bracket
  (only need "a factor divides the product", NOT the spare-factor strengthening — so `Q = β+1` may
  be `< 4`, which is why the subtop does NOT use Lemma F's `Q ≥ 4`).
- `threerow_c2_boundary_top` — top index `Δ(b+2) ≥ 2`. `a` odd: subadditivity + Lemma A (`b+4` odd
  `≥ 7`). `a` even: bridge collapses `Δ = 2[v₂(∏_{t=1}^{β+2}(P+t)) − v₂(P+β+1)]`, Lemma F (`Q ≥ 4`).
- `threerow_c2_boundary_sub` — subtop index `Δ(b+1) ≥ 1`. `a` odd: subadditivity + Lemma A (`b+2`
  odd `≥ 5`). `a` even: the polynomial identity `W₁ = 2(2(β+1)(P+1)+β)` (proved by `ring` after
  substitution) + bridge + `vz_dvd_le`.
- `threerow_c2_boundary` — the assembled statement `V < valb1 ∧ V < valb2` under `V = val(0)`
  (`θ = 0`, since `j₀ = 0` for `c = 2`).

## What remains MN-verified (NOT in Lean, by design)
The two §3 closed forms are taken as the hypotheses `hΔ2`, `hΔ1`:

> `Δ(b+2) = (b+6) + 2 s₂(P) − 2 s₂(a+2) − 2 v₂(a)`,
> `Δ(b+1) = (b+3) + 2 s₂(P+1) − 2 s₂(a+2) + 2 v₂(W₁) − 2 v₂(a)`,  `W₁ = a(b+2) − b(b+1)`.

They rest on the Murnaghan–Nakayama closed forms for `M_0` (hook), `M_{b+2}` (Lemma T) and
`M_{b+1} = (b−1)W₁/2` — symmetric-function facts out of Mathlib's reach. The Lemma-T factor
`v₂(a−b+1) = 0` is folded in (since `a−b+1 = 2P+3` is odd). Same scoping as every prior session.

**Independent numerical check (this session):** the two closed forms were re-verified against the
true MN valuations `val(b+i) − val(0)` (`G_j = C(m,j) M_j`, `val(j) = j + 2 v₂(G_j)`) over `b ∈
[3,21]`, `P ∈ [0,13]` — **532 cases, 0 mismatches**, and the proved bounds (`Δtop ≥ 2`, `Δsub ≥ 1`)
all hold. So the Lean hypotheses are not vacuous.

## Convention note
Per the LEAN brief's warning: the master / direct convention is used — `W₁` (and the implicit
`N_i`) keeps `(a−b+1)`/`(a−c+2)` inside. We parametrise by `P` with `a = 2P + b + 2` (the Lemma-T
regime `P ≥ 0`, i.e. `a ≥ b + 2`). The finitely many `a = b` shapes (`P = −1`) sit outside the
Lemma-T closed-form derivation and are the small-shape casework of the interior note, not the
boundary lemma proper.

## Scope honesty
Formalises the **number-theoretic content only** — the `Δ` reductions as integer/`padicValNat`
arithmetic plus the existing kernels. The `M_j` closed forms and `G_λ` itself (symmetric-function
content) remain out of Mathlib reach, MN-verified rather than Lean-checked, consistent with the
whole project.
