# Lean snapshot — three-row `c = 3` boundary lemma, end-to-end

**Date:** 2026-06-18 (lean session) · **Project:** `lean/tworow_d4_kernel`, Mathlib v4.30.0

`ThreeRowC3Boundary.lean` assembles the boundary lemma for the three-row tie family `λ = (a, b, 3)`
into a single `sorry`-free declaration `threerow_c3_boundary`, closing the boundary residual of the
`c = 3` interior theorem (`2026-06-16-c3-boundary-complete.md`). With `c = 1`
(`2026-06-16-lean-c1-boundary/`) and `c = 2` (`2026-06-18-lean-c2-boundary/`) this makes **all three**
complete three-row families end-to-end machine-checked.

## Status
- `lake build`: 0 errors, 0 warnings (full project: 2091 jobs, success).
- `#print axioms threerow_c3_boundary` → `[propext, Classical.choice, Quot.sound]` only (no `sorry`,
  no extra axioms). Same for the six per-index sub-theorems and the `vz_N2_ge_one` helper.

## What `c = 3` adds over `c = 1, 2` (the genuinely new content)
- **The two-factor engine `lemma_F2` is used in anger.** The `a`-odd cases present *two* simultaneous
  unit deficits — `v₂(a−1)` (from `a−1 = 2(P+β+1)`) and `v₂(P+2)` (from `a−b+1 = 2(P+2)`). The top
  index `b+3` and subtop `b+2` each dissolve both at once with `lemma_F2` (`Q ≥ 6`, surplus `≥ 1`).
  This is the first session where `lemma_F2`'s `Q ≥ 6` threshold is load-bearing.
- **`N₁, N₂` evenness is derived in Lean, not assumed.** From the *definitions*
  `N₂ = a(b+3) − (b²+2b+3)` and `N₁ = ab² + 5ab + 6a − (b³+b²+4b+6)` (taken as hypotheses), the file
  proves `1 ≤ v₂(N₂)` (exact decomposition `N₂ = 2(Pb+3P+2b+3)`, helper `vz_N2_ge_one`) and
  `1 ≤ v₂(N₁)` (parity-split decompositions `N₁ = 2·(…)`, `ring` after substituting `a = 2P+b+3` and
  `b = 2β` / `2β+1`). These plug the deficit the bridge cannot.
- The odd-`c` **convention** matters here (per the LEAN brief): `N_i` is the master/direct form
  keeping `(a−b+1)`/`(a−c+2)` inside. The expanded closed forms below extend the `c = 2`
  hypothesis shape verbatim.

## What is machine-checked here (pure 2-adic number theory)
- `vz_N2_ge_one` — `1 ≤ v₂(N₂)` from the definition of `N₂` and `a = 2P+b+3`.
- `threerow_c3_boundary_top_aeven` / `_aodd` — top index `b+3`: `Δ ≥ 2` (`a` even, `lemma_F`),
  `Δ ≥ −1` (`a` odd, `lemma_F2` on `t = 2, β+1` of `∏_{t=1}^{β+2}(P+t)`).
- `threerow_c3_boundary_sub2_aeven` / `_aodd` — subtop `b+2`: `Δ ≥ 1` (`a` even, `vz_dvd_le` +
  `v₂(N₂) ≥ 1`), `Δ ≥ 0` (`a` odd, `lemma_F2` on `t = 1, β` of `∏_{t=1}^{β+1}(P+1+t)`).
- `threerow_c3_boundary_sub1_aeven` / `_aodd` — subtop `b+1`: `Δ ≥ 2` (`a` even, `lemma_F` gives
  `v₂(∏) ≥ 1`, plus `v₂(N₁) ≥ 1`), `Δ ≥ −1` (`a` odd, `lemma_F` on `t = β−1` of `∏_{t=1}^{β}(P+2+t)`
  plus `v₂(N₁) ≥ 1`).
- `threerow_c3_boundary` — the assembled statement `V < valb1 ∧ V < valb2 ∧ V < valb3` under
  `V = val(0)` (`θ = 0`, `a` even, `b ≥ 7`) or `V = val(0) − 3` (`θ = 3`, `a` odd, `b ≥ 10`).

Reused kernel (from `ThreeRowC2Boundary` / `LemmaF` / `HookKummerLemmas`): `vz_prod_Icc` (the bridge),
`vz_dvd_le`, `lemma_F`, `lemma_F2`, `sum_digits_two_mul`, `sum_digits_two_mul_add_one`.

## What remains MN-verified (NOT in Lean, by design)
The three expanded Lemma D closed forms are taken as the hypotheses `hΔ3, hΔ2, hΔ1`:

> `Δ(b+3) = (b+5) + 2 s₂(P) − 2 s₂(a+2) − 2 v₂(a−1) − 2 v₂(P+2)`,
> `Δ(b+2) = (b+2) + 2 s₂(P+1) − 2 s₂(a+2) + 2 v₂(N₂) − 2 v₂(a−1) − 2 v₂(P+2)`,
> `Δ(b+1) = (b−1) + 2 s₂(P+2) − 2 s₂(a+2) + 2 v₂(N₁) − 2 v₂(a−1)`.

These are the §1 Lemma D forms with the `carries(2(P+k), b+2k±1)` terms expanded via
`carries(x,y) = s₂(x)+s₂(y)−s₂(x+y)` (the `s₂(b+2k±1)` columns cancel). They rest on the MN closed
forms for `M_0` (hook), `M_{b+3}` (Lemma T) and `M_{b+1}, M_{b+2}` (with factors `N₁, N₂`) —
symmetric-function facts out of Mathlib's reach. Same scoping as every prior session. `N₁, N₂`
themselves enter only through their *definitions* (`hN1, hN2`); their valuations are proved.

**Independent numerical check (this session):** the three expanded closed forms were re-verified
against the true MN valuations `val(b+i) − val(0)` (`G_j = C(m,j) M_j`, `val(j) = j + 2 v₂(G_j)`,
via Murnaghan–Nakayama) over the box-interior range `m ∈ [8,30)` — **573 cases, 0 mismatches**. The
six bridge-collapse identities and both `N₁` even-decompositions were verified symbolically (sympy).
So the Lean hypotheses are not vacuous.

## Convention / regime
Parametrised by `P` with `a = 2P + b + 3` (the Lemma-D regime `P ≥ 0`, i.e. `a ≥ b + 3`). The
box-interior thresholds `b ≥ 7` (`a` even) and `b ≥ 10` (`a` odd) keep the 2-adic box inside
`[0, b]`; the finitely many smaller-`b` shapes are the per-family boundary-tie casework of the
interior note, not the boundary lemma proper.

## Scope honesty
Formalises the **number-theoretic content only** — the `Δ` reductions as integer/`padicValNat`
arithmetic, the `N_i` valuations, and the existing kernels. The `M_j` closed forms and `G_λ` itself
(symmetric-function content) remain MN-verified rather than Lean-checked, consistent with the whole
project. No new `sorry`, no gap discovered in the paper proof.
