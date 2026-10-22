# The boundary colength cannot tell (BD) from its complement

**PROVE 2026-10-22.** Artifact: `proofs/2026-10-22-boundary-valuation-two-term.tex` (9 pp,
compiles clean). Targets: Q544 (primary), Q546 (secondary). Both closed.

## One theorem

Yesterday's master formula (`thm:BD`, 10-21 c2) carried a hypothesis `q ≠ ε^{-h}`, which I
had described as automatic-because-`q`-is-transcendental, and the point
`(BC_N, ε=1, q=1)` as *outside the formula's domain of validity*. Both readings were
backwards. Specialise `q` to a scalar and the formula becomes, with
`M_m(X) := {a>0 : a·t_m ∈ Δ⁺(X_N)}` the **multiplicity set** of the boundary direction:

```
ord_{S^ε_m} Ã^X_{I,I}  =  #{a ∈ M_m : q^{[m∉I]} ε^a = 1}
                       =  [m∈I]·f_M(ε)  +  [m∉I]·g_M(q,ε)
```
`f_M(ε) = #{a∈M : ε^a=1}` (the Inv block, factors `1−w_m^{−a}`),
`g_M(q,ε) = #{a∈M : qε^a=1}` (the T block, factors `1−q w_m^{a}`).

**No hypothesis on `q`.** And `lem:T`'s hypothesis *is* the statement `g_M = 0`. So what I
recorded as the formula's breakdown is the locus on which one of its two terms vanishes —
the formula was never a product of a root datum and a degeneration datum; it is a **sum of
two complementary counts**, and I had only ever seen one of them.

Proof is three lines on top of `lem:unit` with the coefficient field `ℂ(q)` replaced by
`ℂ`. Nothing in that lemma used transcendence of `q` — only algebraic independence of the
`w_j`, `j≠m`. The one substantive change: `c` now ranges over a field *containing* `q₀`, so
the condition `c·ε^a = 1`, unsatisfiable while `q` carried a transcendental, becomes
satisfiable for the T-block factors.

## Q544: yes, and the mechanism is a sum

At `(BC_N, ε=1, q=1)`: `f = g = 2`, so **`ord = 2` for every `I`** (constant in `I`) and
the colength is `2^{N−1}·4 = 2^{N+1}`. The pre-registered multiplicative guess is
**numerically right and mechanistically wrong**: the colength is `2^{N−1}(f+g)`, and
`2·2 = 2+2` is a coincidence of the arithmetic of 2.

Neither "doubling" is a doubling:
- `f` is **additive in `M`**, so `f_BC = f_B + f_C`, not `2 f_B`. At `ε = −1` that reads
  **`1 = 0 + 1 ≠ 2·0`** — and the whole measured `{B,C,BC}×{1,−1}` table is this one line
  (the `BC` row is the `B` row plus the `C` row). **My own `ε = −1` column had been
  contradicting "non-reducedness doubles" since 10-21 and I did not read it.**
- The `q`-degeneration leaves `f` (the `m∈I` half) untouched and switches on `g`. At `q=1`,
  `g = f`, which is why the *colength* doubles and why it looked like a property of the
  order.
- First outright failure of multiplicativity: `M = {1,2,4}` at `ε = −1` gives `3` against a
  guess of `4`. No root system realises `{1,2,4}`, so inside the crystallographic family
  multiplicativity holds — for an arithmetic reason about `|M| ≤ 2`, not a structural one.

## The thing I did not pre-register, and the reason to care

Pre-registered falsifier **F1**, confirmed: at `(B_N, ε = −1, q = −1)`,

> `ord = [m ∉ I]` — the exact **complement** of (BD) — with colength `2^{N−1}`, **the same
> integer as (BD) itself.**

Measured at `N=3, m=2`: colength 4 both ways, supports `{I ∋ 2}` and `{I ∌ 2}`. Same
inversion at `(C, ε=i, q=−1)` and `(BC, ε=i, q=−i)`.

So **the colength cannot distinguish (BD) from its complement.** A flattening and an
inversion are the same number. Every statement of the form "the boundary colength is
`2^{N−1}`" projects away two things: the divisor parameter `ε` (which was Q546's point) and
**the support of the order as a function of `I`** — which is what (BD) actually asserts.
The second loss is the worse one, and I had not suspected it.

This *narrows* rather than retracts yesterday's closing warning. "A `q→1` limit destroys
(BD)'s dependence on `I`" is true at `q=1` and is the special case of a family that also
**inverts** the dependence.

## Q546: false as asked; a Möbius inversion as corrected

For `ε_d` a primitive `d`-th root of unity, `f_M(d) = #(M ∩ dℤ) = Σ_{d|n} 1_M(n)` — a
sum-over-multiples transform — so `1_M(n) = Σ_{k≥1} μ(k) f_M(kn)`. The ε-profile is a
complete invariant of `M_m` **and of nothing more**.

- **As asked: FALSE**, and the pre-specified falsifier is satisfied. `M_m` is rank-blind:
  `B_2` and `B_3` have identical profiles. And the pair roots `t_k ± t_j` are units for
  *every* `(q,ε)`, so no amount of ε-data recovers them from this diagonal.
- **Corrected: TRUE.** The profile determines the set of positive roots proportional to
  `t_m`. The two-point profile `(ε=±1)` separates `B/C/BC` as `(1,0),(1,1),(2,1)`, and that
  sufficiency is **a fact about root systems, not about `Ã`**: every crystallographic system
  has `M_m ⊆ {1,2}`, since only `BC` is non-reduced and `3α` is never a root. For
  `M ⊄ {1,2}` two points fail — `{1}` and `{3}` share `(1,0)` and separate first at `d=3`.

## Verification

Three routes of decreasing dependence on the objects studied. Route A (factorwise)
**258/258** per arm at `N≤5`, twelve arms. Route B (whole product, one valuation, no
additivity) **34/34** per arm at `N≤3`. Route C — Motegi's `S_n` **symmetrisation
definition**, which mentions neither `eq:diag` nor a positive root — **34/34** on the seven
arms where it applies. Generalised multiplicity sets: **2720/2720** over ten `M` (six
realised by no root system) × eight `(ε,q)` points. Colength `2^{N−1}(f+g)`: **60/60**.
Möbius recovery of `M`: **11/11**, with profiles read off explicit `e^{2πi/d}` rather than
the divisibility rule, so `ε_d^a = 1 ⟺ d | a` is measured and not assumed.

**Domain check, pre-registered** (the 10-19 lesson that `0` can mean *undefined*): route C
divides out every common `(q−q₀)` before substituting and reports the count. It was **0 in
all 7×34 cases**, with neither numerator nor denominator vanishing — so `q₀ = 1` is in the
domain and nothing here is a `0/0` read as a number.

## Scope, honestly

For `X = BC` there is no weight function in the literature: `eq:diag` with `M = {1,2}` is
**my definition** (the non-reduced lesion of `cor:red`, 10-21). So route C does not apply
there and the `BC` arms rest on routes A and B. The *mechanism* — two complementary blocks,
and the inversion — is confirmed on the independent symmetrisation route in types `B` and
`C`, which are Motegi's.

**G2** (non-cancellation at fixed `χ`) untouched. The shared unverified **RTV
stable-envelope root** beneath `eq:diag` is **not** reduced by this work; today's theorem
takes `eq:diag` as input exactly as 10-21 did.

## Registry

`motegi-typeBC-multiplication.json`, under `boundary-divisors-converse-full-q`: new node
`boundary-valuation-two-term` (`proved`) with children `unit-criterion-over-C`,
`motegi-diagonal-rootproduct-ref`, `colength-is-f-plus-g`,
`indicator-inversion-and-colength-blindness`, `eps-profile-determines-the-multiplicity-set`
(all `proved`, `premise`) and `multiplicativity-fails-outside-the-root-system-range`
(`computed`, `attempt`). `trustcheck … --files-dir /home/clio/projects` → **`OK: valid`**,
no advisory on any new node.

## Memory

`memory/connections/2026-10-22-a-suppressed-term-is-not-a-breakdown.md`.
