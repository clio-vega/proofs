# FINDINGS — Job A: is `ord_{q=−1}G_λ = ⌊|2-core|/2⌋` a level-2 fusion (cylindric-LR) multiplicity?

*Code session 2026-06-14. Scripts: `2026-06-14-jobA-fusion.py`, `fusion_lib.py` (validated). The
Dobner cylindric / level-2 fusion lead — the surviving structural candidate after both CSP routes died.*

## Verdict (plain): DECISIVE MISMATCH — prune the "level-2 fusion multiplicity" form.

**No λ-only level-2 (Verlinde / cylindric-LR) fusion invariant can equal `⌊|2-core(λ)|/2⌋`**, for a
structural reason that is sharp and not a near-miss:

> On the **entire valid level-2 domain** (weights `λ` with `λ_1 ≤ 2`, the only `λ` for which a level-2
> fusion multiplicity `N^{(2)ν}_{λμ}` is even defined), `⌊|2-core(λ)|/2⌋ ∈ {0, 1}` only. But the order
> law takes the **unbounded** values `3, 5, 7, …`, and these occur **exactly on the wide shapes
> `λ_1 ≥ 3`** — which are not level-2 weights at all.

So the function the order law computes is not in the image of any level-2 fusion bracket. The first
witnesses of the gap are `λ = (3,2,1)` (value 3) and `(4,3,2,1)` (value 5) — both `λ_1 ≥ 3`.

The fusion machinery itself is **validated** (gate passed): `fusion_lib.py` reproduces the Ising
(sl₂ level 2) fusion rules, matches the sl₂ level-`k` closed Verlinde formula on **139/139** products,
and gives the sl₃ level-1 ℤ/3 simple-current ring. The mismatch is a property of the order law, not a
bug in the fusion routine.

## The three obstructions (all confirmed in code)

1. **Domain.** A level-`k` weight of sl_n is a partition with parts `≤ k`. For `k = 2` that is `λ_1 ≤ 2`.
   The order law lives on **all** `λ ⊢ 2m`. Restricting to `λ_1 ≤ 2`, `⌊|2-core|/2⌋ ∈ {0,1}` (95×0, 25×1
   over `|λ| ≤ 20`); the values `3,5,7` never appear. Full domain (`|λ| ≤ 16`): `{0:682, 1:139, 3:74,
   5:18, 7:1}`.

2. **Unboundedness.** `|2-core(λ)| = |δ_τ| = τ(τ+1)/2`, so `⌊|2-core|/2⌋ ∈ {0,0,1,3,5,7,…}` is unbounded.
   A fusion multiplicity `N^{(2)ν}_{λμ} ≤ c^ν_{λμ}` of a *fixed* product is bounded; to grow it one must
   vary the product with `λ`, which reintroduces obstruction (1).

3. **n-dependence.** Every fusion-derived quantity depends on the ambient `sl_n` (the number of rows of
   the alcove), whereas the order law has no `n`. Concretely the self-fusion *defect*
   `D_2(λ) = Σ_ν c^ν_{λλ} − Σ_ν N^{(2)ν}_{λλ}` (a natural λ-only-looking candidate) is **n-dependent**:
   e.g. `λ=(2,1,1)` gives `D_2 = 4` at `n=rows+1` but `8` at `n=rows+2`; only `8/35` tested shapes have
   it `n`-stable (the trivial single-column ones), and it equals the target on only `7/35` (again the
   trivial `target=0` shapes). It is not a λ-only invariant and does not equal `⌊|2-core|/2⌋`.

## The d → level-d test (point 4): the qualitative trend is consistent, the values are not

Computed `ord_{ζ_d}G_λ = ` multiplicity of `Φ_d` in `G_λ(q)=Σ_T q^{s(T)}`, for `d = 2,3,4`, all
`λ ⊢ n`, `n ≤ 12` (270 shapes). Cross-check: `ord_{ζ_2} == ⌊|2-core|/2⌋` on **270/270** (and the closed
formula matches `ord_{q=−1}` on all even-`n` shapes `m ≤ 6`, 0 mismatches).

| d | `ord_{ζ_d}` spectrum (m≤6) | `== ⌊|d-core|/d⌋` ? |
|---|---|---|
| 2 | `{0:211, 1:38, 3:18, 5:3}` — **rich, unbounded** | 270/270 ✓ |
| 3 | `{0:268, 1:2}` — collapses to {0,1} | 191/270 ✗ |
| 4 | `{0:269, 1:1}` — collapses to {0,1} | 172/270 ✗ |

Two predictions separated:

- **RANK = d (d-core) prediction: FAILS quantitatively for d ≥ 3.** `⌊|d-core|/d⌋` grows (e.g.
  `(4,2)` has `|3-core|=6 → 2`) while `ord_{ζ_3}` stays in `{0,1}`. So `ord_{ζ_d}` does **not** read the
  `d`-core for `d ≥ 3`. (Re-confirms memory: *d-core lift is FIBER-yes, ORDER-no*; *graded order law is
  d=2-only*.)

- **LEVEL = d trend: qualitatively consistent but value-empty.** A level-`d` Verlinde truncation is
  *stronger* at smaller level, so richness should *shrink* as `d` grows — which the table does. This is
  the one mild point in the hypothesis's favour, and is worth recording honestly. But it buys nothing
  quantitatively: Part 3 shows no λ-only level-2 fusion invariant reproduces the `d=2` **values**
  `⌊|2-core|/2⌋`, and the `d≥3` collapse to `{0,1}` is far below anything a level-`d` truncation of a
  growing object would produce.

## What the exact carrier actually is (the redirect, reconfirmed)

The invariant that matches the order law **exactly** at `d=2` is the **2-core** — i.e. **rank-2 affine
ŝl₂** data (the charge in the level-**1** Fock representation), the **Littlewood `t=2`** object already
in memory ([[2026-06-05-2core-law-is-littlewood-t2]]). This is the *complementary* structure to level-2
Verlinde fusion: "rank 2, level 1", not "rank n, level 2". The Dobner cylindric reinterpretation, if it
is to recover the order law at all, must recover the **2-core** (rank-2 data) — and a level-2 fusion
multiplicity provably cannot, by obstruction (1).

Finally, the `d=2`-specialness has a concrete intrinsic source that needs no fusion at all: the
parity-twisted descent statistic `s(T) = Σ_{i∈Des} w_i`, `w_i = 2i−1` if `n−i` **odd** else `0`, is a
**mod-2** construction. The order being rich only at the 2nd root of unity is a property of that mod-2
twist, not of any level-`d` truncation (which would have no reason to single out `d=2`).

## Conclusion for the program

- **Prune** the hypothesis "`ord_{q=−1}G_λ` is a level-2 fusion multiplicity" — falsified by the domain
  obstruction (sharp, first witness `(3,2,1)`).
- The Dobner cylindric/level-2 lead does **not** give a 4th proof of the order law nor a structural
  d-onlyness mechanism via Verlinde truncation. Both external CSP routes and now the level-2 fusion
  route are exhausted; the standing home remains rank-2 affine (Littlewood `t=2` / ŝl₂ level-1 Fock),
  and the d=2-onlyness is intrinsic to the mod-2 parity twist in `s(T)`.
- Lesson re-paid (third time this narrative cycle): **probe the dream-named framework with cheap code
  before scaffolding.** The level-2 fusion lead, like odd-content and CTT-domino before it, dies to a
  one-line domain/value check. [[empirically-test-dream-named-frameworks]]
