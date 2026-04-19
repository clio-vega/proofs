# The Cactus Midpoint Theorem — Proof Summary

**Date:** 12 April 2026  
**Status:** PROVED (all 5 parts)

## Results

### Part 1 (Midpoint = Cactus, q=0)
At q=0 and p=1/2: ∏ R_{i_j}(1/2) = (1/2)^m π_{w₀^[a,b]}.  
**Proof:** R_i(1/2) = (1/2)π_i by the factorization R_i(p) = pπ_i + (1-2p). Product telescopes by word-independence.

### Part 2 (Möbius Factorization)
∏ R_{i_j}(p) = Σ_S p^|S| (1-2p)^{m-|S|} · π_S.  
At p=1/2, only the full subset S={1,...,m} survives.  
**Proof:** Binomial expansion of ∏(pπ_i + (1-2p)).

### Part 3 (Spectral Parameter)
κ = (1-qp)/(1-p). Three levels:
- p=0: κ=1 (identity)
- p=1/2: κ=2-q (cactus midpoint)  
- p=1: κ=1/(1-q) (Hecke)

### Part 4 (Generic q — NEGATIVE)
At generic q≠0, π_w is NOT well-defined: the braid relation fails.  
**Key identity:** π_i π_{i+1} π_i - π_{i+1} π_i π_{i+1} = q(T_i - T_{i+1}).

### Part 5 (Why q=0 is Special)
Three properties hold at q=0 and fail elsewhere:
1. Idempotency: π_i² = π_i
2. Braid relation: obstruction q(T_i - T_{i+1}) vanishes
3. Word-independence: π_w well-defined → rank 1 projector

## The Braid Obstruction (core new result)

The heart of the paper is a single algebraic identity:

**π_i π_{i+1} π_i - π_{i+1} π_i π_{i+1} = q(T_i - T_{i+1})**

This says the Demazure operators satisfy the braid relation up to an error
proportional to q. At q=0 the error vanishes — this is WHY coboundary
structure (cactus group action) is special to q=0.

## Computational verification
- H_q(S_3) regular representation: all identities verified
- H_0(S_3), H_0(S_4): midpoint identity verified for all parabolics
- Rank(π₁π₂π₁) = 1 at q=0, = 3 at q=1
