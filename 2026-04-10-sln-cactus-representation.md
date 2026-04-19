# sl_n Cactus Representation Theorem — Proof and Correction

**Date:** 2026-04-10  
**Status:** Part 1 PROVED, Parts 2-3 CORRECTED (original claims refuted)

## Theorem (proved)

The reversal permutation σ_{[i,j]} on V^{⊗k} (V = C^n, any n ≥ 2) defines a representation of the cactus group J_k. The proof is dimension-free: it is identical to the sl_2 case.

## Original claims refuted

**Claim 2 (= HK crystal commutor): FALSE.** The Henriques-Kamnitzer crystal commutor on B(ω₁)^{⊗k} with identical factors is the identity (trivial). The reversal permutation is nontrivial and is NOT the same as the HK commutor.

**Claim 3 (intertwines crystal operators): FALSE.** The reversal does NOT commute with Kashiwara crystal operators on a single tensor product crystal structure. Counterexample: for sl_2, f₁(1,1) = (2,1), but swap(f₁(1,1)) = (1,2) ≠ f₁(swap(1,1)) = (2,1).

## Correct crystal-theoretic statement

The reversal σ_{[p,q]} IS a crystal isomorphism between:
- (B(ω₁)^{⊗k}, standard reading order)  
- (B(ω₁)^{⊗k}, reordered reading with positions p..q reversed)

This is tautological: the Kashiwara tensor product rule is defined relative to a reading order, and permuting values together with the reading order preserves the reduced signature.

## Why the HK commutor is trivial for identical factors

When B₁ = B₂ = B(ω₁), the crystal structures on B₁⊗B₂ and B₂⊗B₁ are identical (the tensor product rule depends on φ and ε values, not on copy labels). The decomposition B(ω₁)⊗B(ω₁) = B(2ω₁) ⊕ B(ω₂) is multiplicity-free, so the unique crystal automorphism is the identity.

## Computational verification

- Cactus relations verified for sl_2, sl_3, sl_4 (k=3,4)
- Crystal morphism failure verified for sl_2, sl_3, sl_4
- Crystal intertwining (standard ↔ reordered) verified for all intervals

## Key insight

The reversal representation lives in the category of REPRESENTATIONS (where P is the braiding/commutor), not in the category of CRYSTALS (where the commutor for identical objects is trivial). This is consistent with Alqady-Stroinski: the coboundary monoidal structure lives in TL_0 (a representation category), not in crystals.

## Files

- `2026-04-10-sln-cactus-representation.tex` — full LaTeX proof (compiles to 6 pages)
- `2026-04-10-sln-cactus-representation.pdf` — compiled proof
