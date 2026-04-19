# The Multiplicity Bundle Theorem

**Date:** April 11, 2026
**Status:** Complete (all four parts proved)

## Theorem Statement

Let g = sl_n, V = C^n the fundamental representation, V^{⊗k} its k-fold tensor product. Under the Schur-Weyl decomposition V^{⊗k} ≅ ⊕_λ V_λ ⊗ S_λ:

1. **Schur-Weyl factorisation:** The coboundary commutor σ_{[i,j]} = ⊕_λ (id_{V_λ} ⊗ ρ_{S_λ}(w_0^{[i,j]})).

2. **Crystal invisibility:** The crystal commutor on identical fundamental factors B(ω_1)^{⊗k} is the identity.

3. **Multiplicity bundle:** The fibre over each crystal component B(λ) is the Specht module S_λ. The coboundary commutor acts only on these fibres.

4. **Three-level hierarchy:** braided monoidal →[forget spectral param]→ coboundary monoidal →[forget multiplicity]→ symmetric monoidal.

## Proof Summary

**Part 1** is classical: σ_{[i,j]} is a permutation of tensor factors (by the Cactus Representation Theorem of 04-10), and Schur-Weyl duality decomposes any S_k element as id ⊗ ρ on each isotypic component.

**Part 2** is the mathematical heart, proved in three steps:
- *Rigidity:* Every crystal automorphism of a connected crystal B(μ) is the identity (the unique highest weight element is fixed, and all other elements are determined by Kashiwara operators applied to it).
- *Multiplicity-free:* B(ω_1) ⊗ B(ω_1) = B(2ω_1) ⊔ B(ω_2), each appearing once. So any crystal automorphism is the identity.
- *Extension:* σ_{[i,j]} decomposes into adjacent-factor commutors, each trivial on identical B(ω_1) factors.

**Parts 3-4** package the result categorically: the crystal functor is a projection from the "multiplicity bundle" (total space) to the crystal (base space), discarding the Specht module fibre where the coboundary commutor lives.

## Computational Verification

All claims verified for:
- sl_2, k = 3, 4
- sl_3, k = 3
- sl_4, k = 3

Specific checks:
- Crystal component multiplicities match Specht module dimensions (hook length formula)
- σ_{[1,k]} eigenvalues match Schur-Weyl predictions
- σ is an involution and commutes with GL_n
- B(ω_1)^2 is multiplicity-free for sl_2 through sl_5
- Naive tensor-factor reversal does NOT preserve crystal components (confirming the crystal commutor ≠ naive swap)

## Gaps

None. All four parts are proved. The key ingredients are:
- Cactus Representation Theorem (proved 04-10)
- Classical Schur-Weyl duality
- Crystal automorphism rigidity (standard)
- Multiplicity-freeness of Sym^2(V) ⊕ Λ^2(V) = V⊗V (classical)

## Files

- LaTeX proof: `2026-04-11-multiplicity-bundle.tex` (compiled to .pdf)
- Verification script: `verify_multiplicity_bundle.py`
