# Lean snapshot — c=4 boundary content is axiom-free; durability fixed; redundancy witnesses packaged

**Date:** 2026-07-05 · **Project:** `projects/lean/tworow_d4_kernel` (Mathlib v4.30.0, Lean 4)
· **Target (LEAN.md):** (1) put the Lean sources under version control; (2) make the c=4 boundary
2-adic content lemmas axiom-free (Rick's ask, emails 411/412).

## Headline

- **Durability FIXED.** The live Lean sources — including the previously untracked
  `ThreeRowC4InteriorN4.lean` (the c=4 interior N4 keystone) — are now committed to
  `clio-vega/proofs` under `lean-src/tworow_d4_kernel/` and pushed. No machine-checked work is
  now at risk of being lost.
- **The c=4 boundary content lemmas were ALREADY axiom-free** (derived in the Jun 16–19 cycles,
  not re-discovered this session). Verified this cycle: **no content bound `k ≤ vz N_i` is assumed
  as a hypothesis anywhere in the project** — every one is a *derived theorem* from the explicit
  even decomposition `N_i = 2^{k_i} · (integer)`. `#print axioms` on all of them shows only
  `[propext, Classical.choice, Quot.sound]`.
- **New this session:** two clean *named* packaging theorems giving Rick a single citable
  multiplicative-side counterpart to his additive `additive_redundancy_at_eS` — zero new
  mathematics, just conjoining the existing derived floors.

## What is already derived, axiom-free (confirmed, not assumed)

In `ThreeRowC4Boundary.lean`, each `N_i` content floor is a theorem proved via the
`hNeq : N = 2^k * inner` + `dvd_mul_right` + `padicValNat_dvd_iff_le` pattern (the `⟨_, ring⟩`
witness that FAILS is avoided — explicit `2^k *` form used):

| lemma | statement | decomposition |
|---|---|---|
| `vz_N3_ge_one_beven` | `1 ≤ vz N₃` (b even) | `N₃ = 2^1·(2βP+5β+4P+4)` |
| `vz_N2_ge_two` | `2 ≤ vz N₂` (both parities) | `N₂ = 2^2·(…)`, inner differs by parity |
| `vz_N1_ge_three_beven` | `3 ≤ vz N₁` (b even) | `N₁ = 2^3·(…)` |
| `vz_N1_ge_two_bodd` | `2 ≤ vz N₁` (b odd) | `N₁ = 2^2·(…)` |

These are wired into all eight boundary-index lemmas (`threerow_c4_boundary_{top,sub1,sub2,sub3}_{beven,bodd}`)
and hence into the assembled `threerow_c4_boundary`, which depends on the std three axioms only.

The c=3 boundary (`ThreeRowC3Boundary.lean`) likewise derives its content (`vz_N2_ge_one`, N₁
evenness inline); the c=4 interior `N4` (`16 ∣ H`) is self-contained via `decide` over `ZMod 16`.
So the whole `tworow_d4_kernel` project is **0 sorry, 0 `axiom` declarations, std-three only**.

The remaining *assumed* hypotheses in the assembled theorems (`hΔ4…hΔ1`, and the `hN1…hN3`
polynomial **definitions**) are the Murnaghan–Nakayama closed forms / factor definitions for the
symmetric-function object `G_λ` — genuinely out of Mathlib's reach, and NOT content facts. Rick's
ask was specifically about the content floors `v₂(N_i) ≥ …`; those are theorems.

## New this session: packaged multiplicative-redundancy witnesses

Added to `ThreeRowC4Boundary.lean` (zero new math — `⟨…⟩` over existing derived lemmas):

```lean
theorem multiplicative_redundancy_c4_beven … :
    1 ≤ vz N₃ ∧ 2 ≤ vz N₂ ∧ 3 ≤ vz N₁ := ⟨…⟩
theorem multiplicative_redundancy_c4_bodd … :
    2 ≤ vz N₂ ∧ 2 ≤ vz N₁ := ⟨…⟩
```

These are the multiplicative-side counterpart to Rick's `additive_redundancy_at_eS`: one citable,
axiom-free statement per parity of `b` bundling the fixed 2-adic floors of the redundancy factors,
for the joint FREE/RIGID note.

## Build / axioms

- `lake build`: **green, 2093 jobs, 0 sorry, 0 warnings.**
- `#print axioms`:
  - `threerow_c4_boundary` → `[propext, Classical.choice, Quot.sound]`
  - `N4` → `[propext, Classical.choice, Quot.sound]`
  - `multiplicative_redundancy_c4_beven` / `_bodd` → `[propext, Classical.choice, Quot.sound]`
  - `vz_N2_ge_two`, `vz_N1_ge_three_beven`, `vz_N3_ge_one_beven` → `[propext, Classical.choice, Quot.sound]`

## Commits (pushed to `clio-vega/proofs`)

- `270b348` — add canonical `lean-src/tworow_d4_kernel` tree (durability).
- `5855b2d` — c4 boundary: packaged multiplicative-redundancy content witnesses.

## Honest scoping

- I did **not** invent any new content lemma this cycle because none was missing — the LEAN.md
  belief that "axiom-free content asks were NOT met" was stale for the *content floors themselves*
  (they were met Jun 16–19). The genuinely unmet ask was **durability** (now fixed) and a
  **single named witness** for the joint note (now added).
- No `N_i` is factored; the `c ≥ 4` irreducibility wall is irrelevant to these floors.
