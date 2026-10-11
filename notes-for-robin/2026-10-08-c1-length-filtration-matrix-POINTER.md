# POINTER — PROVE 2026-10-08 c1: what `M^(L)` is

**This directory has no git remote. Robin cannot read it.** The readable note is in the pushed repo:

- https://github.com/clio-vega/proofs/blob/main/2026-10-08-for-robin-length-filtration-matrix.md
- https://github.com/clio-vega/proofs/blob/main/2026-10-08-length-filtration-matrix.tex (9 pp.)

Content: `M^(L)`'s diagonal blocks are *diagonal* (entry `∏_i m_i(ρ)!`), hence
`det M^(L) = ∏_{ℓ(ρ)≤L} ∏_i m_i(ρ)!`, explicit `(M^(3))^{-1}`, Theorem B at `L=3`, and a Bell(k)
bound on how many `c_{λ,ν}` a `k`-part class sees. Plus two verification findings: the determinant is
**provably blind** to the off-diagonal claim it appears to corroborate, and `(-1)**ht` with `ht<0`
returns a Python *float*, which made sympy report 86 true identities as failures.
