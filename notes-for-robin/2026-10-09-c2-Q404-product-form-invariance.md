# For Robin — PROVE 2026-10-09 c2 (Q404): Theorem D now holds for all b ≥ 3, and it is presentation-invariant

**Paper:** `proofs/2026-10-09-c2-product-form-invariance.tex` (13 pp, compiles clean).
**Registry:** `proofs/registry/two-part-green-polynomials.json`, +14 nodes under `thm-D-product-form-obstruction`.

## The one-line mathematical content

A root of `D_b = t^b − t^{b−1} + 1` on the unit circle must satisfy `|α − 1| = 1` *as well as*
`|α| = 1`, because `D_b(α) = 0` says exactly `α^{b−1}(α − 1) = −1`. Two unit circles meet in two
points. So at most two roots of `D_b` lie on the circle, and since `deg D_b = b > 2` for `b ≥ 3`,
some root is off it.

That replaces the real-root argument, which was **provably** unavailable for even `b` — for even `b`,
`D_b` has no real roots at all. The published Theorem D's hypothesis "`b` odd ≥ 3" can be deleted.

Two by-products of the same lemma, which I did not have before:

- **`Φ₆` is the only cyclotomic polynomial that can divide `D_b`**, and it does exactly when
  `b ≡ 2 (mod 6)`, with multiplicity one (every root of `D_b` is simple). This is the sharp form of
  the warning not to claim irreducibility: `D_8 = Φ₆·(t⁶−t⁴−t³+t+1)`.
- The companion family `E_b = t^b − t^{b−1} + 2` (the `a = b` diagonal) has **no root in the open unit
  disc**, because `|α| < 1` gives `|α^{b−1}(α−1)| < 2`. Hence `M(E_b) = ν(E_b) = 2` *exactly*, for
  every `b ≥ 2`.

## What I think is the better result

The whole two-row × two-part slice is now **completely classified**. With `λ = (a,b)`, `ρ` two-part,
`m = min(ρ)`:

> `Y^{(a,b)}_ρ` admits an admissible product form **iff** `m ≠ b`, or `b = 1`, or (`b = 2` and `a ≥ 3`).

Every off-diagonal cell *is* a product form, explicitly:

- `m > b`:  `Y = −t^{b−1}(1−t)`
- `1 ≤ m ≤ b−1`:  `Y = −t^{b−m−1}(1−t)·(1−t^{2m})/(1−t^m)`

I had never written those down, and they localise the obstruction: it is caused entirely by the
`⟦y=b⟧ + ⟦x=b⟧` term of Theorem C — the only place that formula contributes a *constant*. The three
exceptional product forms are `t`, `Φ₆`, `Φ₂`.

A pleasant tidy-up: the fourth of the four original non-cyclotomic witnesses, `t³+t²−1`, which the
10-07 scope corollary recorded as "occurs only at `ℓ(λ)=3`", is just `−D₃(−t)`. Same invariant
`ν = θ₀² = 1.754877666`. It was never a separate phenomenon.

## Q404 — the answer, and the reason

**Yes, the obstruction is presentation-invariant**, and the reason is structural rather than lucky.

1. *Product form* means membership in the **subgroup** `Ĉ = ℚ^× · ⟨t, 1−t^d (d ≥ 1)⟩` of `ℂ(t)^×`.
   Membership in a subgroup is a property of the **element**, so the predicate cannot be changed by
   rewriting the target. Theorem D was already stated as an existential over expressions, hence as
   non-membership in that set — so it was function-level all along.
2. The obstruction is **certified** by a group homomorphism
   `ν(f) = ∏_{α ≠ 0} max(|α|, |α|^{−1})`, which vanishes on `Ĉ`. `ν` is additionally invariant under
   `t → t^k`, `f → f*`, `k`-th roots, and `t → −t`. So **your widening of the citation to signed
   exponents, products and quotients costs nothing** — that is exactly (R2) in the paper, and it is
   free because `Ĉ` is a group.
3. The class is **saturated**: if `g ∈ ℤ[t]` is monic and not in `±C`, then `ν(g) > 1`. So **no monic
   integer polynomial can be adjoined to the admissible class without introducing a root off the unit
   circle.** `ker ν` is the largest class for which the certificate works.

The disanalogy with the signedness ceiling that died this cycle is now precise. Pak–Robichaux's
laundering `f = [g + (2^{Cn^a} − h)] − 2^{Cn^a}` destroys an inference from a *presentation* ("my
formula has minus signs") to a *function-level* claim. It has no analogue here because `#P` is a
**semigroup** inside the group it generates, while `Ĉ` **is** a group with an invariant on the
quotient. The criterion I take away:

> *An obstruction "`f` admits no expression of shape `S`" is presentation-invariant iff the set of
> values of shape-`S` expressions is specified without reference to `f`; and it is certifiable iff
> that set lies in the kernel of a homomorphism out of the ambient structure that does not kill `f`.
> Group structure on the admissible class is what turns the first into the second.*

## Two corrections to my own work, for the record

- **The invariant I proposed for myself was wrong by one leading coefficient.** I had written the
  target restatement as "no product form ⟺ `M > 1`" with `M` the Mahler measure. That is false once
  rational scalars are admitted: `M(2t) = 2 > 1` while `2t` is a monomial. `M` sees the leading
  coefficient; `ν` does not. On monic targets they agree (`ν = M²` when `|f(0)| = 1`) and `D_b`, `E_b`
  are monic, so **Theorem D is unaffected** — but the invariant correct for the class *as stated* is
  `ν`, and that is what the paper uses.
- **A symmetric test cannot see a reversal.** My first Green-polynomial engine computed *cocharge*
  rather than charge (`K_{λλ}` read `t^{n(λ)}` and `K_{(n),μ}` read `1` — exactly swapped). After I
  corrected the index rule, the *asymmetric* test `deg K_{λμ} = n(μ) − n(λ)` with leading coefficient
  `1` still failed on **133 of 471** nonzero pairs. All three *symmetric* tests I had written
  (`K_{λλ}=1`, `K_{(n),μ}=t^{n(μ)}`, `K(1)=#SSYT`, 434 pairs) read **0 failures** throughout, because
  each is invariant under the convention I had wrong. I abandoned charge and rebuilt the engine from
  the definition (Gram–Schmidt for `⟨p_λ, p_μ⟩_t = δ z_λ ∏(1−t^{λ_i})^{−1}`, triangular over
  dominance). Nothing in the paper depends on charge.

## Verification summary

- Independent engine (definition-based, no shared code path with the 10-07 derivation):
  `Y^μ_ρ(0) = χ^μ_ρ` against Murnaghan–Nakayama, **434 pairs `n ≤ 7`, 0 failures**, control `Y→Y+1`
  fires 434/434; `Y^μ_ρ(1) = [m_μ]p_ρ` by direct count, **434 pairs, 0 failures**; Theorem C
  reproduced on all **28** two-row × two-part cells `n ≤ 7`, **0 failures**, with a control dropping
  `⟦a=b⟧` predicted in advance to fire on 3 cells and firing on exactly that set.
- Exact (non-numerical) cyclotomicity test, self-tested in both arms: `D_b` is a monomial times
  cyclotomics precisely for `b ≤ 2` and never for `3 ≤ b ≤ 40`; `E_b` precisely for `b = 1`; all
  **465** off-diagonal cells `b ≤ 30` pass.
- `ν(D_b) = M(D_b)²` for `2 ≤ b ≤ 60`, 0 failures. `M` computed two ways (roots, and Jensen's contour
  integral) agreeing to `1.6e−6`, the residual being quadrature error — the first comparison reported
  "36 of 40 failures" and that was the *instrument*, not the value (log singularity on the contour
  exactly when `Φ₆ | D_b`).
- `pdflatex` 0 errors, 0 undefined references; planted-error control reads 1 and 1.
- `trustcheck` 0 problems; two-arm control on the same file fires (broken path → 1, bogus enum → 2).

## One thing I did not do, and one honest gap

- **I did not email Rick.** His FPSAC deadline is 2026-11-15 and commit `d7bca5e5b` widened
  `[Clio26, Thm D]` to signed exponents, products and quotients — which §(R2) now shows is free, and
  the `b ≥ 3` improvement strengthens what he is citing. This was a prove session (no email by the
  session rules), so the note is here and the paper is pushed; sending it is the obvious first action
  for the next session that may.
- **Smyth (1971) is cited from memory and was not read.** It would give the uniform bound
  `M(D_b) ≥ θ₀ = 1.3247…` for all `b ≥ 3`, since I prove the hypothesis `D_b* ≠ ±D_b`
  (`D_b* − D_b = t^{b−1} − t`). The measurement matches to **43 digits** at `b = 3`. The node is
  graded `speculative` and **no theorem in the paper depends on it.** Likewise the Lawton-type limit
  `M(D_b) → exp(m(1+x+y)) ≈ 1.38136` is recorded as an observation only (measured ≈ 1.3815 for
  `b ∈ [50,60]`).
