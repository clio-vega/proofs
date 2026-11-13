# Q585 H1 — Purbhoo's Cor 3.4 swap M vs tableau switching Σ (PROVE 2026-11-13) — NOT PROVED; reduced

**Theorem attempted.** M = Σ on LR(ν/μ; λ) (β_CM inputs). Status: computed 894/894 distinct elements (376 new this
session, census max(|λ|,|μ|)=6, c≥2); not proved.

**Proved (from Purbhoo Thm 2 / Lemma 2.5, verified-quote).**
1. Factorization: M = Φ_G ∘ κ^{-1}, κ = Cor 3.2 bijection (Table 1 row 1), Φ_G = "flock (−E_B,−N_B) on B, migrate B→A".
   Hence H1 ⟺ Φ_G(P) = Σ(κ(P)) for every mosaic P — two readouts of ONE puzzle.
2. Rotation covariance: Φ_G(P) = R_iii(ρ^{-1}P) (definitional).
3. c=1 cells: trivial.

**Refuted.** The brief's local lemma (one migration step = one jdt move of Σ): G-journey k pairs with slide |μ|+1−k
(reverse order, 74/74), and the hexagon step-string is not a function of the slide move-string ('R': 20 step-strings).
The local lemma must be stated on a puzzle readout — that is Purbhoo Prop 5.2.

**Reduction (gap).** Given Prop 5.1 (computed 134/134) and K1: κ = hw∘revrot∘Tao (computed 134/134, 67/67 cells),
H1 ⟺ Tao∘ρ^{-1} = Σ∘K∘Tao — the form of Cor 5.3 ("By Proposition 5.2"). Remaining: (G1) Prop 5.1, (G2) K1,
(G3) Cor 5.3/Prop 5.2 — Purbhoo's omitted proofs; (G4) pure jdt identity f↦g̃ = (Σ∘K)^{±1} — next attackable step.
Tie-break: still an interpretation (12431 ties, 0 unresolved).

**H2:** not-located — reading blocked by session rule (no browsing); 0906.0077 abstract states reducibility, not equality.

Artifact: `2026-11-13-Q585-H1-migration-is-switching.tex/.pdf`; code `code-1113-Q585-H1/` @ 29b0d30e9c4c.
