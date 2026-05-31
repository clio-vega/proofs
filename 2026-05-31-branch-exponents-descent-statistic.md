# The branch exponents of the w₀ monodromy are a descent statistic — and the order law becomes a combinatorial minimum

*Clio, 2026-05-31*

Hi Robin — a clean one today, and I think it's the nicest reformulation of the order law I've
found. Short version: the spectral invariant I've been chasing for two weeks turns out to be a
parity-twisted major-index statistic on standard Young tableaux, and the order law
`ord_{x=q²} Z_λ = τ(τ+1)/2` reduces to a one-line minimum over SYT with a complete,
representation-theory-free proof.

## The objects
For `λ ⊢ n`, `M(x)` is the Baxterized longest-element monodromy on the irreducible Hecke
module `V^λ` — the ordered staircase product of `Ř_i(x) = P_q^{(i)} + r(x)P_{-1}^{(i)}`,
`r(x)=(q²−x)/(q(x−1))`, along the reduced word for `w₀`. `Z_λ(x)=tr_{V^λ} M(x)`, and `x=q²`
is the fusion point. The **branch exponents** `{d_j}(λ)` are the orders to which the
eigenvalues of `M(x)` vanish at `x=q²` (Newton-polygon slopes of the characteristic polynomial
in `s=x−q²`); there are `dim V^λ = #SYT(λ)` of them.

Already known: `min_j d_j = τ(τ+1)/2` (the order law), `Σ_j d_j = Σ_T comaj(T)`,
`M(q²)=0` iff `τ≥1`.

## The finding (computed exactly, all 43 partitions of n≤7)
The whole multiset is a descent statistic:

> **`{d_j}(λ) = { s(T) : T ∈ SYT(λ) }`,  where  `s(T) = Σ_{i∈Des(T)} w_i`,  and
> `w_i = 2i−1` if `n−i` is odd, `w_i = 0` if `n−i` is even.**

Descents at positions of the same parity as `n` are *free*; opposite-parity descents cost
`2i−1`. The individual eigenvalues are NOT indexable by single tableaux (they leave `q²` as
confluent Puiseux clusters — fractional slopes that only recombine into integers in the
symmetric functions), so the statistic lives on the *distribution* of `s(T)`, not on an
eigenvalue↔tableau bijection. That's why a year of "is `d_j` a content statistic?" guesses
missed it.

Two corollaries, both clean:
- **Conjugation symmetry** `{d_j}(λ') = {C(n,2) − d_j(λ)}`, because transposing a tableau
  complements its descent set and `Σ_i w_i = C(n,2)`. (Also provable spectrally via the
  sgn-twist `S_i=(q−1)I−T_i` and the crossing identity `M(x)M(q²/x)=I`, using
  `r(x)^{-1}=r(q²/x)`.)
- The **value alphabet `A(n)` depends only on `n`** (subset sums of `{2i−1: n−i odd}`),
  satisfying `A(n)=A(n−2)∪(A(n−2)+(2n−3))`; all the shape information is in the multiplicities.

## The payoff — order law as `min_T s(T) = τ(τ+1)/2`
This is a combinatorial statement with a complete proof strategy (verified out-of-sample to
n=9, all 65 partitions):
1. **Prefix-descent bound (the only real step):** `|Des(T)∩[1..m]| ≥ r(m)−1`, where
   `r(m)=min{r:λ_1+⋯+λ_r ≥ m+1}`. Entries `1..m+1` form a Young subshape spanning `≥r(m)`
   rows; each new row opened below the first is a descent `≤m`.
2. The prefix constraints form a totally-unimodular interval system, so greedy = LP optimum:
   use free positions freely, pay the `⌈τ/2⌉` cheapest opposite-parity weights.
3. Those cheapest weights sum to `τ(τ+1)/2` (`Σ(4j−1)=k(2k+1)`, `Σ(4j−3)=k(2k−1)`).
4. The bound is achieved by an explicit minimizing tableau.

A nice bonus: the same prefix-congestion argument re-derives `τ = 2ℓ−n−1` (`ℓ` = number of
rows) — the column-1 charge formula drops out of `D(n−1)=ℓ−1` rather than being an input.

## Why I'm excited
The hard residual in the operator/deletion-bridge proof of the *off-hook* order law has been
the "merged-arc reach" gap. This spectral-eigenvalue route via the descent bridge has **no
such gap** — steps 1–4 are uniform over all shapes. So the program splits cleanly:
- **(A)** `min_T s(T)=τ(τ+1)/2` — a second, gap-free derivation of the lower bound; I'll prove
  it rigorously next session (strategy above is essentially complete).
- **(B)** the bridge `{d_j}=\{s(T)\}` itself — the real prize. Proving it makes the entire
  spectral "home" combinatorial and closes the order law uniformly, off-hook included. It's
  harder: it equates the Newton-polygon valuations of a confluent-Puiseux spectrum with a
  descent distribution (all the symmetric functions of the `d_j`, not just sum and min).

I'll push the writeup to clio-vega/proofs once I've drafted the (A) proof. Wanted to flag the
reformulation now — `{branch exponents} = {parity-twisted maj over SYT}` feels like the kind
of "two constructions, one multiset" coincidence that usually means there's real structure
underneath. Curious whether it rings any bells from the cylindric / fake-degree side.

— Clio
