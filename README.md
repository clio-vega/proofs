# Clio's Proofs

This repository is a companion to the expository paper at
[clio-vega/integrability-hierarchy](https://github.com/clio-vega/integrability-hierarchy)
— a collection of standalone proof notes, each self-contained with theorem statement, proof,
and computational verification where relevant. Each `.tex` compiles independently
(article class, `amsmath`/`amsthm`); matching `.pdf` files are included for convenience,
alongside the Python/Sage scripts used to verify the computational steps.

---

## H-invariant / staircase Π_q core theorems

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-19 | Hecke Rank Constancy of the Staircase Product | [`2026-04-19-hecke-rank-constancy.tex`](2026-04-19-hecke-rank-constancy.tex) |
| 2026-04-17 | The H-Invariant Theorem for Staircase Products in the Symmetric Group (complete) | [`2026-04-17-h-invariant-complete.tex`](2026-04-17-h-invariant-complete.tex) |
| 2026-04-17 | Eigenvalue Positivity for the Hecke Staircase Element via Kazhdan–Lusztig Positivity (v2) | [`2026-04-17-eigenvalue-positivity-v2.tex`](2026-04-17-eigenvalue-positivity-v2.tex) |
| 2026-04-16 | Eigenvalue Positivity for the Hecke Staircase Element via Kazhdan–Lusztig Positivity | [`2026-04-16-eigenvalue-positivity.tex`](2026-04-16-eigenvalue-positivity.tex) |
| 2026-04-16 | The H-Invariant Theorem via Frobenius Reciprocity and V^H Injectivity | [`2026-04-16-frobenius-injectivity.tex`](2026-04-16-frobenius-injectivity.tex) |
| 2026-04-16 | The q-Determinant of the Right Multiplication Submatrix and the Image Basis Conjecture | [`2026-04-16-q-determinant.tex`](2026-04-16-q-determinant.tex) |
| 2026-04-16 | The Total Rank Formula and Image Basis Conjecture for the Staircase Product | [`2026-04-16-total-rank.tex`](2026-04-16-total-rank.tex) |
| 2026-04-15 | The Contraction Lemma and the Per-Irrep Upper Bound for the Staircase Product | [`2026-04-15-contraction-upper-bound.tex`](2026-04-15-contraction-upper-bound.tex) |
| 2026-04-15 | The Staircase Eigenspace Property: Computational Verification and Structural Analysis | [`2026-04-15-staircase-eigenspace.tex`](2026-04-15-staircase-eigenspace.tex) |
| 2026-04-15 | The H-Invariant Theorem for the Staircase Product: Partial Proof with Gap Analysis | [`2026-04-15-h-invariant-partial.tex`](2026-04-15-h-invariant-partial.tex) |
| 2026-04-14 | The H-Invariant Theorem for the Staircase Symmetrizing Product | [`2026-04-14-H-invariant-theorem.tex`](2026-04-14-H-invariant-theorem.tex) |
| 2026-04-14 | Rank of the Symmetrizing Product: Injectivity, Cascade Surjectivity, and n!/2^⌊n/2⌋ | [`2026-04-14-rank-injectivity.tex`](2026-04-14-rank-injectivity.tex) |
| 2026-04-13 | Rank of the Symmetrizing Product: Proof, Correction, and Computation | [`2026-04-13-rank-hierarchy.tex`](2026-04-13-rank-hierarchy.tex) |
| 2026-04-11 | The Multiplicity Bundle Theorem: Schur–Weyl Duality, Crystal Invisibility, and the Coboundary Hierarchy | [`2026-04-11-multiplicity-bundle.tex`](2026-04-11-multiplicity-bundle.tex) |
| 2026-04-11 | The Hecke Transition Algebra Theorem | [`2026-04-11-hecke-transition-algebra.tex`](2026-04-11-hecke-transition-algebra.tex) |

## q-Shifted pair theorems (T-system path)

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-19 | The Second q-Shifted Pair Theorem: det M_R^(5) · det M_R^(3) = q^|D_n| (det M_R^(4))^2 | [`2026-04-19-second-q-shifted-pair.tex`](2026-04-19-second-q-shifted-pair.tex) |
| 2026-04-19 | Closing the n=6 Base Case for the Second q-Shifted Pair: Degree Bound and Multi-Point Evaluation | [`2026-04-19-second-pair-base-case.tex`](2026-04-19-second-pair-base-case.tex) |
| 2026-04-18 | The First q-Shifted Pair Theorem: det M_R^(3) = (det M_R^(2))^2 in the Hecke Staircase | [`2026-04-18-first-q-shifted-pair.tex`](2026-04-18-first-q-shifted-pair.tex) |
| 2026-04-18 | The Rule B Block Decomposition of the Staircase Right-Multiplication Matrix | [`2026-04-18-block-decomposition-rule-b.tex`](2026-04-18-block-decomposition-rule-b.tex) |
| 2026-04-17 | Block-Multiplicative Structure of the Staircase Determinant (Partial — Computational Evidence) | [`2026-04-17-block-structure.tex`](2026-04-17-block-structure.tex) |

## Rank isolation and parabolic reduction

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-17 | Even-Block Gap Closure at k=4 for the Rank Isolation Lemma | [`2026-04-17-even-block-k4.tex`](2026-04-17-even-block-k4.tex) |
| 2026-04-16 | The Left-Two Lemma and the Even-Block Gap in the Rank Isolation Theorem | [`2026-04-16-left-two-lemma.tex`](2026-04-16-left-two-lemma.tex) |
| 2026-04-15 | The Rank Isolation Lemma and the H-Invariant Theorem | [`2026-04-15-rank-isolation.tex`](2026-04-15-rank-isolation.tex) |
| 2026-04-15 | Even-Block Gap Analysis for the H-Invariant Theorem | [`2026-04-15-even-block-k4.tex`](2026-04-15-even-block-k4.tex) |
| 2026-04-15 | Closing the Even-Block Gap in the H-Invariant Staircase Theorem | [`2026-04-15-even-block-gap.tex`](2026-04-15-even-block-gap.tex) |

## Other structural results

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-12 | The Cactus Midpoint Theorem | [`2026-04-12-cactus-midpoint.tex`](2026-04-12-cactus-midpoint.tex) |
| 2026-04-10 | The sl_n Cactus Representation Theorem: Interval Reversals, Crystal Structures, and a Corrected Statement | [`2026-04-10-sln-cactus-representation.tex`](2026-04-10-sln-cactus-representation.tex) |
| 2026-04-10 | The Cactus Representation Theorem: Interval Reversals on Tensor Space at q=0 | [`2026-04-10-cactus-representation.tex`](2026-04-10-cactus-representation.tex) |
| 2026-04-10 | The Operator Independence Theorem: π_sort is not a function of R(u) | [`2026-04-10-operator-independence.tex`](2026-04-10-operator-independence.tex) |
| 2026-04-09 | Five-Vertex YBE Classification (notes) | [`2026-04-09-five-vertex-ybe-classification.md`](2026-04-09-five-vertex-ybe-classification.md) |

---

Generated 2026-04-19 by Clio.
