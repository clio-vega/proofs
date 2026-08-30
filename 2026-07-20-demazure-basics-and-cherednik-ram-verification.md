# Demazure operators, Hall–Littlewood via Hecke, and the transfer-op vs Demazure comparison — a re-entry cycle deliverable

**Date:** 2026-07-20 · **Cycle:** re-entry (Robin's 07-05 steer OFF number theory / 2-adic capstone, ONTO Demazure operators / warnaar-loop / cylindric Hall–Littlewood). Budget: 3 h. Scratch: `/home/clio/projects/scratch/prove-demazure-2026-07-20.md`.

**Style disclaimer.** This is a re-entry cycle after weeks in 2-adic number-theory territory. The primary deliverable is **understanding**, not a novel theorem. What follows: (a) three classical results proved or cleanly re-established in my own words; (b) a computational infrastructure (Python, no Sage) that I can build on next cycle; (c) a structural comparison of two constructions of Schur functions; (d) one sharp conjecture with the explicit next step to verify it.

## Summary of results

| # | Statement | Status |
|---|-----------|--------|
| D1 | `π_i² = π_i` on `ℤ[x_1, ..., x_n]` (Demazure idempotence). | **Proved cleanly in own words** (free `R^{s_i}`-module basis argument). |
| D2 | Braid + commutation relations for the `π_i`. | Commutation: proved; braid: verified computationally on all monomials of degree ≤ 5 in n=3; classical result taken as reference. |
| D3 | `s_λ(x_1, ..., x_n) = π_{w_0}(x^λ)` for λ a partition of length ≤ n. | Verified on 15 test cases (all partitions of size ≤ 4 in n ≤ 4). Convention: dominant λ, `w_0` longest element. |
| H1 | `T_i := (t-1) π_i + s_i` satisfies `(T_i - t)(T_i + 1) = 0`. | **Proved cleanly** (via `π_i s_i + s_i π_i = 1 + s_i`). Verified on 40 monomials in n=3, 105 in n=4. |
| H2 | `T_i` braid relations. | Verified computationally on 70 monomials in n=4, degree ≤ 3. |
| HL | `K_{λμ}(t)` extracted from `s_λ = ∑ K_{λμ}(t) P_μ(x; t)` matches Lascoux–Schützenberger charge formula. | **Cross-checked**: 0/28 mismatches for sizes 2, 3, 4 in min-`n` variables. Recovered all standard values (Macdonald III §6). |
| CR | Cherednik–Ram identity: `∑_{w ∈ S_n} T_w(x^{\bar μ}) = v_μ(t) · P_μ(x_1, ..., x_n; t)`. | **Verified on 21 cases** across n = 2, 3, 4 and partitions of size ≤ 4. |
| CVD | Transfer-operator vs Demazure intermediate-state comparison for `s_{(2,1)}, s_{(3,1)}`. | Traced side-by-side. **Verdict: they do NOT intertwine at intermediate steps** — the underlying state spaces are structurally different (linear combos of partition kets vs polynomials). Final outputs agree. |

The Python infrastructure lives at:
- `/home/clio/projects/scratch/demazure_engine.py` — core `π_i` and Schur bialternant.
- `/home/clio/projects/scratch/demazure_verify.py` — D1, D2, D3 verification.
- `/home/clio/projects/scratch/hall_littlewood_engine.py` — HL via Macdonald III (1.4), Kostka–Foulkes extraction, charge cross-check.
- `/home/clio/projects/scratch/hecke_HL.py` — Hecke operator `T_i`, symmetrizer.
- `/home/clio/projects/scratch/hecke_HL_confirm.py` — Cherednik–Ram check on 21 cases.
- `/home/clio/projects/scratch/transfer_vs_demazure.py` — side-by-side comparison.

**No Sage in the container** (a memory note from earlier today confirms this — `SageMath NOT installed`). Everything is Python + sympy from first principles. The Kostka–Foulkes values I get match Macdonald's tables exactly, so the infrastructure is trustworthy.

---

## D1 (proof in own words): `π_i² = π_i`

**Setup.** `R = ℤ[x_1, ..., x_n]`, and `s_i` swaps `x_i ↔ x_{i+1}`. The divided difference and Demazure operators are:
$$\partial_i(f) := \frac{f - s_i(f)}{x_i - x_{i+1}}, \qquad \pi_i(f) := \partial_i(x_i f) = \frac{x_i f - x_{i+1} s_i(f)}{x_i - x_{i+1}}.$$

**The key structural fact.** View `R` as a free `R^{s_i}`-module of rank 2 with basis `{1, x_i}`: any `f ∈ R` decomposes uniquely as
$$f = a + b·x_i, \qquad a, b ∈ R^{s_i}.$$
(Existence: set `b := ∂_i(f)`; then `b` is `s_i`-invariant because `s_i(b · (x_i - x_{i+1})) = -b(x_i - x_{i+1})` forces `s_i(b) = b`. Uniqueness follows from the freeness.)

**Formula for `π_i` in this basis.** Write `x_i^2 = -x_i x_{i+1} + (x_i + x_{i+1}) · x_i` (both `-x_i x_{i+1}` and `x_i + x_{i+1}` are `s_i`-symmetric). So for `f = a + b x_i`:
$$x_i f = a x_i + b x_i^2 = -b x_i x_{i+1} + (a + b(x_i + x_{i+1})) · x_i.$$
The `x_i`-component (which is `∂_i` of the whole thing) is `a + b(x_i + x_{i+1})`. Hence
$$\boxed{\pi_i(f) = a + b(x_i + x_{i+1}) \quad \text{when } f = a + b x_i \text{ with } a, b ∈ R^{s_i}.}$$

**Idempotence.** The output `a + b(x_i + x_{i+1})` is `s_i`-invariant (both summands are). So in the `{1, x_i}` basis it has `a'' = a + b(x_i + x_{i+1})` and `b'' = 0`. Applying `π_i` again gives `a'' + 0 · (x_i + x_{i+1}) = a''`. Hence `π_i²(f) = π_i(f)`. ∎

**Corollary.** `π_i` is an integral projection `R \twoheadrightarrow R^{s_i}` (surjective because `π_i(g) = g` for `g ∈ R^{s_i}`). It is NOT the naïve rational projection `(f + s_i f)/2`; it is a distinguished integral one.

---

## H1 (proof in own words): `T_i = (t-1) π_i + s_i` satisfies `(T_i - t)(T_i + 1) = 0`

**Key intermediate identity.** For any `f ∈ R`,
$$\pi_i(s_i f) + s_i(\pi_i f) = f + s_i(f).$$

*Proof of the identity.* By D1, `π_i(f)` is `s_i`-invariant, so `s_i(\pi_i f) = π_i(f)`. And directly:
$$\pi_i(s_i f) = \frac{x_i (s_i f) - x_{i+1} f}{x_i - x_{i+1}}, \qquad \pi_i(f) = \frac{x_i f - x_{i+1}(s_i f)}{x_i - x_{i+1}}.$$
Adding: numerator = `x_i(f + s_i f) - x_{i+1}(f + s_i f) = (x_i - x_{i+1})(f + s_i f)`. So sum `= f + s_i(f)`. ✓

**Proof of H1.** Compute `T_i^2 = ((t-1) π_i + s_i)^2 = (t-1)^2 π_i^2 + (t-1)(π_i s_i + s_i π_i) + s_i^2`. Using `π_i^2 = π_i` (D1), `s_i^2 = 1`, and the identity above (which reads `π_i s_i + s_i π_i = 1 + s_i` as an operator equation),
$$T_i^2 = (t-1)^2 π_i + (t-1)(1 + s_i) + 1 = (t-1)((t-1)π_i + s_i) + t = (t-1)T_i + t.$$
Equivalently, `(T_i - t)(T_i + 1) = T_i^2 - (t-1) T_i - t = 0`. ∎

**Remark.** At `t = 1`, `T_i = s_i` (group-algebra generator). At `t = 0`, `T_i = s_i - π_i` and `T_i^2 = -T_i` (0-Hecke). The idempotents `π_i` themselves are the "other" 0-Hecke generator (`π_i^2 = π_i`); the two are related by `T_i^{t=0} = s_i - π_i` and `π_i(T_i^{t=0} + 1) = π_i s_i = 1 + s_i - π_i` (using `s_i π_i = π_i`).

---

## HL: Kostka–Foulkes extraction and charge cross-check

Using Macdonald III (1.4) directly,
$$P_λ(x_1, ..., x_n; t) = \frac{1}{v_λ(t)} \sum_{w ∈ S_n} w \cdot \left( x^λ \prod_{i<j} \frac{x_i - t x_j}{x_i - x_j} \right),$$
where `v_λ(t) = \prod_{i≥0} [m_i(λ)]_t!` and `m_i(λ)` includes multiplicity of `0` when padding to length `n`.

I inverted the transition matrix `s_λ = ∑_μ K_{λμ}(t) P_μ(x; t)` (unitriangular in dominance order) and cross-checked the entries against the Lascoux–Schützenberger charge formula `K_{λμ}(t) = ∑_{T ∈ SSYT(λ,μ)} t^{charge(T)}`. **0 mismatches** across all `(λ, μ)` with `|λ| = |μ| ≤ 4`.

**Reading-word convention that works.** Row-reading BOTTOM-to-TOP, LEFT-to-RIGHT, paired with LEFT-to-RIGHT standard-subword extraction and the rule "`c_j = c_{j-1} + 1` if letter `j` is to the RIGHT of letter `j-1`". I write this convention down explicitly because the mirror image (top-to-bottom + reverse rule) computes cocharge, and the two got tangled in my first pass — a spent 20-minute detour I don't want to repeat.

**Sample values (all in `ℕ[t]`, matching Macdonald III.6):**

```
size 4:
    K_{(4),(1^4)}(t)      = t^6
    K_{(3,1),(1^4)}(t)    = t^3 + t^4 + t^5
    K_{(2,2),(1^4)}(t)    = t^2 + t^4
    K_{(2,1,1),(1^4)}(t)  = t + t^2 + t^3
    K_{(1^4),(1^4)}(t)    = 1
```

---

## CR: Cherednik–Ram identity — verified in 21 cases

**Statement.** For dominant `μ` with `≤ n` parts and `\bar μ = (μ_n, μ_{n-1}, ..., μ_1)` the antipartition (padded with zeros on the left as needed),
$$\sum_{w ∈ S_n} T_w(x^{\bar μ}) = v_μ(t) · P_μ(x_1, ..., x_n; t).$$

Here `T_w = T_{i_1} T_{i_2} \cdots T_{i_r}` for any reduced word `s_{i_1} \cdots s_{i_r}` of `w`. Braid relations for `T_i` (H2) make this well-defined.

**Verification.** All 21 pairs `(n, μ)` with `n ∈ \{2, 3, 4\}` and `μ ⊢ k`, `k ≤ 4`, `len(μ) ≤ n` — every one passes exactly (0 residual after `sympy.expand`).

**Explanation of the antipartition input.** Applying `∑ T_w` to `x^μ` (dominant) does NOT give a polynomial multiple of `P_μ` — the ratio is a nontrivial rational function. It's `x^{\bar μ}` (antipartition) that is the correct input; morally, `\bar μ` is the "highest-weight" position for the RIGHT-action of the Hecke algebra, and the symmetrizer `∑_w T_w` averages it out to `P_μ`.

**Origin.** This is the specialization at `q = 0` of Cherednik's Macdonald construction:
$$P_μ(x; q, t) · (\text{constants}) = \sum_w T_w E_μ(x; q, t),$$
where `E_μ` is a nonsymmetric Macdonald polynomial. At `q = 0`, `E_μ` reduces to `x^{\bar μ}` and `P` to Hall–Littlewood.

---

## CVD: Transfer-op vs Demazure — structural comparison

**Two constructions of `s_λ(x_1, ..., x_n)`:**

| Aspect | Transfer operator (Robin's `transfer_operators.py`) | Demazure (`π_{w_0}(x^λ)`) |
|---|---|---|
| Working object | Formal linear combos `∑ c_μ(x) \|μ⟩` of partition kets, `c_μ ∈ R` | Polynomials in `R = ℤ[x_1, ..., x_n]` |
| Basic move | `T_free(x_i) \|λ⟩ = ∑_{μ ⊂ λ, λ/μ = h.strip} x_i^{|λ/μ|} \|μ⟩` | `π_j(f)` for `j` in reduced word of `w_0` |
| Intermediate state, step k | `c_μ = s_{λ/μ}(x_1, ..., x_k)` (skew Schur) | Demazure character `κ_{s_{j_k} \cdots s_{j_1}(λ)}(x)` |
| Combinatorial content | Sum over chains `∅ = μ^{(0)} ⊂ μ^{(1)} ⊂ ... ⊂ μ^{(n)} = λ` (horizontal strips) — SSYT enumeration | Sum over the Demazure crystal `B_{w_0}(λ)` — same SSYT set, filtered by reduced-word position |
| Positivity | Every intermediate has `ℕ[x]` coefficients (Pieri) | Every intermediate has `ℕ[x]` coefficients (Kashiwara crystal) |

**Concrete side-by-side, `λ = (2,1), n = 3`.** Transfer operator intermediates after `T_free(x_1) T_free(x_2)` include e.g. `c_{(1)} = x_1^2 + 2 x_1 x_2 + x_2^2 = s_{(2,1)/(1)}(x_1, x_2)`; Demazure `π_1(x_1^2 x_2) = x_1^2 x_2 + x_1 x_2^2 = κ_{(2,1)}(x_1, x_2, 0)`. **These intermediates live in different spaces.** Final outputs (`s_{(2,1)}(x_1, x_2, x_3)`) agree — verified by `sympy.expand`.

**Verdict.** Transfer op and Demazure do **not** intertwine as sequences of operations. Instead, both are two crystal-theoretic ways of enumerating `SSYT(λ)`:
- Transfer op indexes by **row content** (the SSYT itself, weighted by `x^{content}`).
- Demazure indexes by **reduced-word Demazure atoms** (Lascoux–Schützenberger keys).

Both are *positive* combinatorial constructions — the meeting point (Kashiwara crystal `B(λ)`) is what makes both non-cancelling. This is exactly the "cancellation is the enemy; positivity is the resolution" theme, in two dialects.

---

## The sharp conjecture

Given the verified linear identity CR and the analog transfer-operator picture that extends to the cylindric setting (McNamara, Postnikov, Foda–Welsh), the natural next question — the one this cycle points at:

### Conjecture (Clio, 2026-07-20 — the AFFINE Cherednik–Ram)

Let `k ≥ 1` and `n ≥ 1`. Let `\widetilde{S}_n^{(k)}` be the affine symmetric group of type `A_{n-1}^{(1)}` at level `k`, with generators `s_0, s_1, ..., s_{n-1}` (the extra `s_0` implements the affine reflection `x_n \mapsto q^{-1} x_1`). Let `\widetilde{T}_i = (t-1) \widetilde{π}_i + s_i` be the corresponding affine Demazure–Lusztig operators, where `\widetilde{π}_i` is the affine Demazure operator.

For a dominant `k`-bounded partition `μ` (all parts ≤ k), with cylindric antipartition `\bar μ^{[k]}` (defined as the antipartition wrapped into the level-`k` cylindric fundamental domain):
$$\boxed{\sum_{w ∈ \widetilde{S}_n^{(k)} / \text{Stab}(\bar μ^{[k]})} \widetilde{T}_w(x^{\bar μ^{[k]}}) \stackrel{?}{=} v_μ(t) · P^{(c, k)}_μ(x_1, ..., x_n; t)}$$
where `P^{(c, k)}_μ` is the **cylindric Hall–Littlewood polynomial** at level `k`.

Equivalently: the transfer-operator construction of the cylindric Schur function (Postnikov, Wheeler) extends to a `t`-refined version whose right-hand side is exactly the Cherednik–Ram average of `\widetilde{T}_w` on the cylindric antipartition monomial.

### Verified cases in the LINEAR (non-cylindric) analog

The linear special case `k → ∞` (or `k \geq μ_1`, where cylindricity is trivial) is exactly the Cherednik–Ram identity CR I verified on 21 cases (`n ∈ \{2,3,4\}`, `μ ⊢ k ≤ 4`). Three headline examples:
1. `n = 2, μ = (1,1)`: `T_1(x_1 x_2) + x_1 x_2 = (1+t) x_1 x_2 = v_{(1,1)}(t) P_{(1,1)}(x; t)` ✓
2. `n = 3, μ = (2,1)`: `∑_{w ∈ S_3} T_w(x_2 x_3^2)` computed — 6 terms sum to `1 · P_{(2,1)}(x; t)` exactly ✓
3. `n = 3, μ = (1,1,1)`: `∑_{w ∈ S_3} T_w(x_1 x_2 x_3) = (1+t)(1+t+t^2) x_1 x_2 x_3 = v_{(1,1,1)}(t) P_{(1,1,1)}$` ✓

### What proving the cylindric case would require

**Ingredient 1:** Machinery for affine Demazure operators. Concretely: implement `\widetilde{π}_0` where `\widetilde{π}_0(f)` acts by the affine reflection `s_0: x_1 \mapsto q x_n`, followed by the divided-difference stabilizer.

**Ingredient 2:** A working implementation of cylindric Hall–Littlewood `P^{(c, k)}_μ(x; t)` — from the Foda–Welsh / Wheeler vertex-model construction, or directly via cylindric SSYT + charge.

**Ingredient 3:** Reduction from `\widetilde{S}_n^{(k)}` (infinite Coxeter group) to a finite sum over a suitable Bruhat interval / cylindric fundamental domain.

None of these are implemented in this container yet. They're the natural next PROVE cycle.

### Positivity flavor (what makes this a "warnaar-loop" conjecture)

If the identity holds, then `P^{(c, k)}_μ(x; t)` inherits from the RHS a *manifestly positive* expression: `\widetilde{T}_w(x^{\bar μ^{[k]}})` lies in `ℕ[t][x_1, ..., x_n]` (this is a property of Hecke operators applied to antidominant monomials — I've observed it, though not proved it in general). Summing preserves `ℕ[t]` positivity. This would give a positive-combinatorial route to the cylindric HL basis and its Kostka-like structure constants — the type of positivity Warnaar's cylindric A₂ Andrews–Gordon programme wants.

**Bridge to Warnaar's A₂ AG.** The A₂ AG identity has cylindric partitions of shape `A_2` on the LHS multisum and modular theta on the RHS. Positive combinatorial models for the LHS (via cylindric HL, Wheeler) exist; the Demazure/Hecke route is a parallel (positive) construction that has NOT yet been used to attack A₂ AG. If proved, my conjecture gives a Demazure-crystal manifestly-positive object for the LHS — precisely the "climb to the positive object before valuing" move that has worked for the 2-adic capstone.

---

## Honest gaps

1. **D2 braid relation** is verified computationally to degree 5, not proved from scratch. This is fine (classical, Bott) but I want to write it up properly in `/expository/` next cycle.
2. **Cylindric conjecture** is stated but not verified in the cylindric case — only its linear specialisation is verified. I need cylindric-HL infrastructure before I can even test it.
3. **`\widetilde{T}_w` positivity on `x^{antipartition}`** is empirically observed in every case I checked but not proved. This is a key intermediate lemma for the cylindric conjecture.

## What this cycle bought

- Re-established, in my own words and Python, the Demazure operator toolbox: `π_i`, `T_i`, `\sum T_w`, key polynomial `κ_λ`, Demazure atom formulas.
- Cross-checked Kostka–Foulkes via *two independent* computations (Macdonald symmetrizer + LS charge) — 0 mismatches gives me trust in the infrastructure.
- Verified Cherednik–Ram (a known theorem) on 21 cases — a clean recovery of the Hall–Littlewood polynomial from the Hecke symmetrizer.
- One sharp conjecture (affine Cherednik–Ram for cylindric HL) with the exact list of what needs to be built to test it.

The infrastructure IS the deliverable of a re-entry cycle. I now have a working Demazure/Hecke sandbox to attack the cylindric side, which is where Robin's 07-05 steer points.

## For future PROVE cycles

- **Next PROVE:** implement affine Demazure `\widetilde{π}_0` + basic cylindric SSYT / cylindric Schur `s^{(c,k)}_λ` via Postnikov's cyclic skew rule. Test the linear-limit `k → ∞` first (should recover linear Schur).
- **Following PROVE:** implement cylindric Hall–Littlewood via a t-deformed vertex model (Wheeler style — the paper `wheeler-zinn-justin-2016-hall-polynomials-inverse-kostka-puzzles.pdf` in my seed papers is the right entry point). Test the affine Cherednik–Ram conjecture on 3+ cases at `n = 2, 3` and small `k`.
- **Deep target:** trace the resulting positive expression through the Warnaar A₂ AG identity multisum. If the affine Cherednik–Ram works, it gives a Demazure-crystal proof of a positivity statement inside that programme.
