# For Robin — PROVE 2026-10-08 c1: what `M^(L)` is

**Readable artifacts (my `memory/` disk has no git remote, so this file lives in the
pushed repo on purpose):**

- https://github.com/clio-vega/proofs/blob/main/2026-10-08-length-filtration-matrix.tex
- https://github.com/clio-vega/proofs/blob/main/2026-10-08-length-filtration-matrix.pdf (9 pp.)
- Registry: https://github.com/clio-vega/proofs/blob/main/registry/two-part-green-polynomials.json
  (13 new nodes under `what-is-M-L`)

## The one-sentence version

Yesterday I proved that the Gram matrix `M_{ρν} = ⟨p_ρ, h_ν⟩` on partitions of `n` is lower
**block**-triangular for length, so every leading block `M^(L)` is invertible. Today's question was
what `M^(L)` *is*. The answer is that the filtration is far finer than block-triangular:

> **each diagonal block is a *diagonal* matrix**, with entry `∏_i m_i(ρ)!` at `ρ`.

The whole proof is one sentence — *a surjection between finite sets of equal cardinality is a
bijection* — and everything below falls out of it.

## What that buys

| | |
|---|---|
| `det M^(L)` | `∏_{ℓ(ρ)≤L} ∏_i m_i(ρ)! = ∏_{ℓ(ρ)≤L} z_ρ/(ρ_1⋯ρ_ℓ)` |
| `det M^(L) = 1` | iff `L=1`, or `n=1`, or (`L=2` and `n` odd) |
| `(M^(L))^{-1}` | `= (M^{-1})^(L)` — the filtration survives inversion |
| row denominators | row `ν` of `M^{-1}` has l.c.m. of denominators **exactly** `∏_k m_k(ν)!` |
| `(M^(3))^{-1}` | explicit; the `ℓ(ρ)=3 × ℓ(ν)=2` block is `−(m_{ν₁}(ρ)+m_{ν₂}(ρ))/(d_ρ(1+⟦ν₁=ν₂⟧))` |
| Theorem B at `L=3` | `Y^λ_{(x,y,z)} = c_{λ,(n)} + Σ_i (1+⟦2ρ_i=n⟧) c_{λ,sort(ρ_i,n−ρ_i)} + (∏_k m_k(ρ)!) c_{λ,ρ}` |
| general shape | support of row `ρ` = **coarsenings** of `ρ`, so a `k`-part class sees at most **Bell(k)** of the `p(n)` coefficients `c_{λ,ν}`. `B₂=2` is Theorem B; `B₃=5`. |

The denominator line is the one I care about most. It says the length-filtration duality is **free in
one direction and costs exactly the automorphism orders in the other** — no primes beyond those
dividing `∏_k m_k(ν)!` ever appear — and by the unimodularity criterion it is an isomorphism of
*integer* lattices only in three degenerate cases.

## Honesty about novelty

The `m`-to-`p` transition (Möbius over the partition lattice) is classical, as is `⟨p_ρ,h_ν⟩` as an
ordered-set-partition count and the determinant of the `S_n` character table — and by `M(p,s) =
M(p,m)M(m,s)` with `|det M(m,s)| = 1`, that last one is *equivalent* to the `L=n` case of my
determinant formula. **I claim no novelty for those three.** No browsing this session, so that is a
recollection, not a literature check. What is new is the *filtered* picture: diagonal-not-just-
invertible, every truncation's determinant, the exact denominators, `(M^(3))^{-1}`, Theorem B at
`L=3`, and the Bell bound.

## Two things in the verification worth your eye

**1. The determinant is blind to the claim it looks like it corroborates.** I checked `det M^(L)`
against exact Gaussian elimination for every `(n,L)` with `n ≤ 9`, zero mismatches, and that reads
like strong corroboration of the structure theorem. It is not. Once triangularity holds, `det` is a
function of the **diagonal alone** — I confirmed this by planting a nonzero in a same-length
off-diagonal slot and watching `det` stay at `24`. A single perturbation of a diagonal block leaves
it triangular. So the determinant agreement is evidence for the *diagonal entries* only, conditional
on the off-diagonal vanishing; the only evidence for the off-diagonal vanishing is the 218 direct
entry checks. I needed a **symmetric pair** of plants to move the determinant at all (`24 → 12`).
This is one rung past "a pass count does not report the rank of the test": here the statistic is
*provably* blind to the sub-claim, by a theorem about determinants.

**2. An instrument fault that read as a refutation.** The first run of the `L=3` Theorem B check
reported **86 of 115 failures**. The polynomials were equal; the comparison was not. In my
Murnaghan–Nakayama recursion the sign was `(-1)**ht` with `ht` a *negative* displacement, and in
Python `(-1)**(-1)` is the **float** `-1.0`, so `sympy`'s `!=` separated `t³+t²+t-1.0` from
`t³+t²+t-1`. Every value correct, every comparison wrong. Had I read that as data I would have
recorded a true theorem as refuted at `n=4`.

## What I did *not* do

- No closed form for `c_{λ,ν}` at `ℓ(ν)=3`. Theorem B at `L=3` reduces the three-part Green
  polynomial to five `h`-coefficients; it does not evaluate them. The missing input is a three-point
  analogue of Rick's two-point shuffle `Sh_{A,B}`, and whether it composes to a triple shuffle is
  untouched.
- Nothing formalised in Lean.
- No literature check (prove session: no browsing).
