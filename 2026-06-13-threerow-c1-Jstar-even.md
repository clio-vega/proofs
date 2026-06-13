# Three-row shapes `(a,b,1)`: `J* ⊆ {j₀, j₀+2}`, so `|J*|` is even on every tie (d = 4)

**Date:** 2026-06-13 (prove session)
**Status:** Sub-family `c=1` (three rows, last row a single box) **PROVED** for all interior `j`
and for the two boundary-tie families; one same-parity boundary inequality verified `m ≤ 80`
(flagged §7). General three-row: exact closed form for `M_j` derived; the precise wall located.

Companion to `2026-06-12-hook-Jstar-even.md` and `2026-06-12-tworow-Jstar-even.md`. Same setup,
same conclusion shape (`|J*|` even on ties), but now the closed form for `M_j` is a **signed
multi-binomial** rather than a single dimension — the first place the two-row `s₂`-collapse
genuinely breaks, and we identify exactly what replaces it.

---

## 0. The problem

Fix `m ≥ 1` and `λ ⊢ 2m`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`,

> `M_j = ⟨s_λ , h₁^{2(m−j)} e₂^j⟩ ∈ ℤ_{≥0}`,  `M₀ = f^λ`,
> `val(j) = j + 2 v₂( C(m,j) M_j )` (`= ∞` when `C(m,j)M_j = 0`),
> `V = min_j val(j)`,  `J* = { j : val(j) = V }`.

Leading-π dichotomy (proved earlier, Lean-checked): with `π = 1+i`,
`G_λ(i) = Σ_j C(m,j) M_j (i−1)^j` and `C ≡ |J*| (mod π)`, so **`|J*|` odd ⟹ `G_λ(i) ≠ 0`**; the
residual ties are exactly the shapes with `|J*|` even.

Throughout `s₂(n)` is the binary digit sum, `v₂` the 2-adic valuation, `v₂(n!) = n − s₂(n)`
(Legendre), `v₂ C(N,k) = s₂(k)+s₂(N−k)−s₂(N)` (Kummer), and `C(n,k) = 0` for `k<0` or `k>n`.

**Main result.** For every three-row shape `λ = (a,b,1) ⊢ 2m` (`a ≥ b ≥ 1`, `a+b = 2m−1`):

> **Theorem.** Exactly one of `a,b` is even, and
> - **`a` even** (so `b` odd): `J* ⊆ {0,2}`, and the tie `J* = {0,2}` occurs **iff
>   `a ≡ 0` and `b ≡ 1 (mod 4)`**;
> - **`a` odd** (so `b` even): `J* ⊆ {3,5}`, and the tie `J* = {3,5}` occurs **iff
>   `a ≡ 3` and `b ≡ 0 (mod 4)`**.
>
> In all cases `|J*| ≤ 2`, so **`|J*|` is even on every tie** (the leading-π layer cancels).

So the offset `j₀ ∈ {0,3}` is governed by the parity of `a`. The interior of the argument
(`0 ≤ j ≤ b−1`) is proved unconditionally below; the two boundary-tie families (`b=1` for `a`
even, `b=4` for `a` odd) are proved in §6; one residual boundary non-tie inequality is verified
`m ≤ 80` (§7).

---

## 1. The general three-row closed form

With `c = 1+x`, expand `(h₁² + x e₂)^m = Σ_j C(m,j) x^j e₂^j h₁^{2(m−j)}`, so for any `λ`
`Σ_j C(m,j) M_j x^j = ⟨s_λ, (h₁²+xe₂)^m⟩`, and `C(m,j) M_j = ⟨s_λ, e₂^j h₁^{2(m−j)}⟩`.

Work in `N = #rows = 3` alphabet variables `s,t,u`, where `h₁ = e₁ = s+t+u` and
`e₂ = st+su+tu`. For a partition `μ = (p,q,r)`, `⟨F, h_p h_q h_r⟩ = [s^p t^q u^r] F(s,t,u)`
(coefficient of the monomial), and `F` is symmetric in `s,t,u`. Define the trivariate coefficient

> `D_j(p,q,r) := [s^p t^q u^r]\, e₁^{2(m−j)} e₂^j
>            = Σ_{i+k+l=j} \dfrac{j!}{i!\,k!\,l!}\, \binom{2(m−j)}{\,p−i−k,\ q−i−l,\ r−k−l\,}`

(a sum of trinomial coefficients; symmetric in `p,q,r`). The three-row Jacobi–Trudi determinant
`s_{(a,b,c)} = det(h_{λ_i−i+j})` gives, with `ℓ = (a−1,\,b−2,\,c−3)`,

> **Proposition 1 (general three-row `M_j`).**
> `  M_j = Σ_{σ∈S₃} sgn(σ)\, D_j\big(a−1+σ(1),\ b−2+σ(2),\ c−3+σ(3)\big).`

This is the **signed multi-binomial**: a `3×3` determinant of trinomial-coefficient sums. The
two-row collapse (single binomial `f^{(a−j,b−j)}`) used the fact that in two variables `e₂ = st`
is a single monomial and `S₂` has two terms; here `e₂` has three monomials and `S₃` has six terms,
so neither the inner `(i,k,l)`-sum nor the outer alternation telescopes. *(Verified against the
symmetric-function definition via Murnaghan–Nakayama for all three-row `λ`, `m ≤ 8`.)*

---

## 2. The `c = 1` closed form

For `λ = (a,b,1)` (`a+b = 2m−1`) the third Jacobi–Trudi index `c−3 = −2` kills four of the six
`S₃` terms (those needing a negative `u`-exponent), and the surviving `D_j` have `r ∈ {0,1}`, which
collapse to ordinary binomials. Writing `N = 2(m−j)`, `β = b−j`:

> **Lemma 2 (`c=1` closed form).** For all `j`,
> `  M_j = b\,C(N,β+1) − (b−j)\,C(N,β) − (j−1)\,C(N,β−1).`
> Two further equivalent forms hold:
> 1. `M_j = C(N,β) + b\,f^{(a−j,\,b−j+1)} + (j−1)\,f^{(a−j+1,\,b−j)}` (with `f^{(p,q)}=C(p+q,q)−C(p+q,q−1)`);
> 2. (for `1 ≤ j ≤ b−1`, i.e. `β ≥ 1`) the **factored form**
>    `  M_j = C(N,\,b−j−1)\,(a−b+1)\,\dfrac{\,b(a+1) − j(j−1)\,}{(b−j)(b−j+1)}.`
> Boundary: `M_b = 2b(m−b)`, `M_{b+1} = b`, `M_j = 0` for `j ≥ b+2`.

*Proof.* The `S₃` reduction gives `M_j = D_j(a,b,1) − D_j(a+1,b−1,1) − D_j(a,b+1,0) +
D_j(a+2,b−1,0)`. For `r = 0`, only the `(i,k,l)=(j,0,0)` term of `e₂^j` survives, giving
`D_j(p,q,0) = C(N, p−j)`. For `r = 1`, the three contributing `(i,k,l)` patterns give
`D_j(p,q,1) = j\,C(N,p−j) + j\,C(N,p−j+1) + N\,C(N−1,p−j)`. Substituting and using `N C(N−1,A) =
(A+1)C(N,A+1)` collapses the `j`-linear pieces and yields the first displayed form; the standard
relations `(N−β)C(N,β) = (β+1)C(N,β+1)` and `(N−β+1)C(N,β−1)=β C(N,β)` give the other two. The
factored form is the polynomial identity `β(β+1) M_j = C(N,β−1)\,B` with
`B = b\,α(α+1) − (b−j)(α+1)(β+1) − (j−1)β(β+1)`, `α = a−j+1` (note `α+β = N`); `B` **factors**:

> **`B = (a−b+1)\,\big(b(a+1) − j(j−1)\big).`**   *(symbolic expansion — §1 of `bexpand.py`)*

Since `α + β = N`, `α = N−β = a−j+1`, the binomial `C(N,β−1)=C(N,b−j−1)`. ∎

*Verified:* all three equivalent forms and the factorization against the symmetric-function
definition, all `λ=(a,b,1)`, `m ≤ 12`.

The factor `(a−b+1)` is **independent of `j`**, so it contributes equally to every `val(j)` and
**cancels in `val(j) − val(j')`.** The factor `b(a+1) − j(j−1)` is the new object: a **quadratic in
`j` inside a valuation**, not a falling-factorial ratio. This is exactly where two rows had nothing
and three rows have the whole difficulty.

---

## 3. The Prop-2 analogue (a clean `val` formula)

Set `K := b(a+1)`. From the factored form, for `1 ≤ j ≤ b−1` (so `M_j > 0`),
`v₂(M_j) = v₂C(N,b−j−1) + v₂(a−b+1) + v₂(K − j(j−1)) − v₂\big((b−j)(b−j+1)\big)`.
A Legendre/Kummer reduction (the falling-factorial ratio `C(2m,b−1)/C(2(m−j),b−j−1) =
(2m)_{2j}/[(b−1)_j(a+2)_j]`, with `v₂((2m)_{2j}) = j + v₂((m)_j)`) collapses everything to:

> **Proposition 2 (`c=1`).** For `2 ≤ j ≤ b−1`,
> `  Δ(j) := val(j) − val(0)
>        = (j−4) − 2 s₂(j−2)
>          + 2\big[\, v₂C(a+2,\,j) + v₂C(b−1,\,j−2) + v₂(b+1) + v₂\big(K − j(j−1)\big) − v₂(a+1)\,\big].`
> Separately `Δ(1) = −1 + 2 v₂(a+2) + 2 v₂(b+1)` and
> `Δ(2) = −4 + 2 v₂(a+2) + 2 v₂(b+1) + 2 v₂(K−2)`.

*(Verified against direct `val(j)−val(0)`, all `λ=(a,b,1)`, `m ≤ 15`, 0 violations.)*

Compare the two-row Prop-2, `val(j)−val(0) = j + 2[v₂C(b,j)+v₂C(a+1,j) − s₂(j)]`. Two-row had a
single bare `s₂(j)`; three-row (`c=1`) has, in addition, the two terms
`+ v₂(K−j(j−1)) − v₂(a+1)` — **the quadratic-valuation term and its `a+1` counterweight.** Since
`a+b = 2m−1` is odd, exactly one of `a,b` is even, and:

- **`a` even** ⟹ `a+1` odd ⟹ `v₂(a+1)=0` and `K = b(a+1)` odd ⟹ `K−j(j−1)` odd ⟹
  `v₂(K−j(j−1)) = 0`. **Both new terms vanish.** (This is the `{0,2}` family.)
- **`a` odd** ⟹ `a+1` even, `b` even ⟹ `K ≡ 0 (mod 4)`; the new terms are active. (The `{3,5}`
  family.)

---

## 4. Two lemmas (one new, two recalled)

> **Lemma A (digit-sum envelope, recalled from the two-row note).** `2 s₂(k) ≤ k+1` with equality
> iff `k ∈ {1,3}`; `2 s₂(k) = k` iff `k ∈ {0,2}`; `2 s₂(k) ≤ k−1` for all `k ≥ 4`.

> **Lemma B (small binomials mod 2, recalled).** `C(n,2)` is odd iff `n ≡ 2,3 (mod 4)`; `C(n,3)`
> is odd iff `n ≡ 3 (mod 4)`.

> **Lemma C (the compensation lemma — NEW, this is the crux).** For all `a,b` and all `1 ≤ j`,
> `  v₂C(a+2,\,j) + v₂\big(b(a+1) − j(j−1)\big) \ \ge\ \ v₂(a+1).`

*Proof of C.* Let `A = v₂(a+1)`. First, the standard inequality `v₂C(n,k) ≥ v₂(n) − v₂(k)`
(from `k\,C(n,k) = n\,C(n−1,k−1)`). Since `a+2` is odd, `j\,C(a+2,j) = (a+2)\,C(a+1,j−1)` gives
`v₂C(a+2,j) = v₂C(a+1,j−1) − v₂(j) ≥ (A − v₂(j−1)) − v₂(j) = A − v₂(j(j−1))`.
Now two cases on `e := v₂(j(j−1))`:
- if `e ≥ A`: then `b(a+1)` and `j(j−1)` are both divisible by `2^A` (as `v₂(b(a+1)) ≥ A`), so
  `v₂(b(a+1)−j(j−1)) ≥ A`, and `v₂C(a+2,j) ≥ 0`; sum `≥ A`. ✓
- if `e < A`: then `v₂(b(a+1)−j(j−1)) = e` (strict-min rule), and `v₂C(a+2,j) ≥ A − e`; sum
  `≥ (A−e) + e = A`. ✓ ∎

*(Lemma C verified for all `λ=(a,b,1)`, `1 ≤ j ≤ b−1`, `m ≤ 30`; the standard inequality for all
`n,k ≤ 200`.)*

Write `R(j) := v₂C(a+2,j) + v₂(K−j(j−1)) − v₂(a+1)`. **Lemma C says `R(j) ≥ 0`.** With it,
Proposition 2 reads, for `2 ≤ j ≤ b−1`,

> `  Δ(j) = (j−4) − 2 s₂(j−2) + 2\big[\, v₂(b+1) + v₂C(b−1,j−2) + R(j)\,\big],`   `R(j),\,v₂C(b−1,j−2) ≥ 0`.

---

## 5. Proof of the Theorem — interior `0 ≤ j ≤ b−1`

Recall `a+b = 2m−1` odd, so exactly one of `a,b` is even. We bound `Δ(j) = val(j)−val(0)` and use
`v₂(b+1) ≥ 1 ⟺ b` odd `⟺ a` even.

### 5.1 Case `a` even (`b` odd) — claim `J* ⊆ {0,2}`, `j₀ = 0`

Here `v₂(b+1) ≥ 1`. We show `Δ(j) > 0` for `j ∉ {0,2}` and `Δ(2) ≥ 0`.

- **`j = 1`:** `Δ(1) = −1 + 2v₂(a+2) + 2v₂(b+1) ≥ −1 + 2 + 2 = 3 > 0` (`a` even ⟹ `v₂(a+2) ≥ 1`).
- **`j = 2`:** `K` odd ⟹ `v₂(K−2) = 0`, so `Δ(2) = −4 + 2v₂(a+2) + 2v₂(b+1) ≥ 0`, with `= 0` iff
  `v₂(a+2)=v₂(b+1)=1`, i.e. **`a ≡ 0` and `b ≡ 1 (mod 4)`** (Lemma B for `C(·,2)`).
- **`j = 3`** (`k=1`): `(j−4)−2s₂(j−2) = −3`; `v₂C(b−1,1)+v₂(b+1) = v₂(b−1)+v₂(b+1) ≥ 3` (two
  consecutive even numbers, one `≡0 (mod 4)`). So `Δ(3) ≥ −3 + 2·3 = 3 > 0`.
- **`j = 4`** (`k=2`): `(j−4)−2s₂(2) = −2`; `C(b−1,2)=(b−1)(b−2)/2` with `b−1` even, `b−2` odd, so
  `v₂C(b−1,2)+v₂(b+1) = v₂(b−1)−1+v₂(b+1) ≥ 3−1 = 2`. So `Δ(4) ≥ −2 + 2·2 = 2 > 0`.
- **`j = 5`** (`k=3`): `(j−4)−2s₂(3) = −3`; `b−1` even ⟹ `b−1 ≢ 3 (mod 4)` ⟹ `C(b−1,3)` even
  (Lemma B), so `v₂C(b−1,3)+v₂(b+1) ≥ 1+1 = 2`. So `Δ(5) ≥ −3 + 2·2 = 1 > 0`.
- **`j ≥ 6`** (`k = j−2 ≥ 4`): `2s₂(k) ≤ k−1` (Lemma A), so `(j−4)−2s₂(j−2) = (k−2)−2s₂(k) ≥ −1`;
  with `2v₂(b+1) ≥ 2`, `Δ(j) ≥ 1 > 0`.

Hence on `0 ≤ j ≤ b−1`, `val` is minimized only at `j ∈ {0,2}`, with `j=2` tying iff `a ≡ 0,
b ≡ 1 (mod 4)`.

### 5.2 Case `a` odd (`b` even) — claim `J* ⊆ {3,5}`, `j₀ = 3`

Here `v₂(b+1) = 0` (`b` even). Plugging `a` odd into the `j=1,2,3` formulas (using `v₂(a+2)=0`,
`v₂(b+1)=0`, `v₂(K−2)=1`, `v₂(K−6)=1`):

`Δ(1) = −1`,  `Δ(2) = −4 + 2 = −2`,  `Δ(3) = (3−4)−2s₂(1) + 2[v₂C(a+2,3) + 0 + 0 + v₂(K−6) − v₂(a+1)]`.
Now `a+2` odd with `a+1` the only even neighbour gives `v₂C(a+2,3) = v₂(a+1) − 1` (Lemma-B
computation), and `v₂(K−6) = 1`, so `Δ(3) = −3 + 2[(v₂(a+1)−1) + 1 − v₂(a+1)] = −3`. Thus

> `val(1) = val(0)−1,\quad val(2) = val(0)−2,\quad val(3) = val(0)−3,`

a **forced descent**: `j = 0,1,2` are never minimal, and the baseline is `j = 3`. Put
`\tilde Δ(j) := val(j) − val(3) = Δ(j) + 3`. For `j ≥ 4`, using `v₂(b+1)=0` and Lemma C
(`R(j) ≥ 0`):

> `  \tilde Δ(j) = (j−1) − 2 s₂(j−2) + 2\big[v₂C(b−1,j−2) + R(j)\big] \ \ge\ (j−1) − 2 s₂(j−2) = (k+1) − 2 s₂(k),`  `k = j−2 ≥ 2`.

By Lemma A: `k=2 (j=4)`: `3−2 = 1 > 0`; `k=3 (j=5)`: `4−4 = 0`; `k ≥ 4 (j ≥ 6)`: `2s₂(k) ≤ k−1`
gives `≥ 2 > 0`. Hence on `0 ≤ j ≤ b−1` the only minimizers are `j = 3` and possibly `j = 5`, with

> `  \tilde Δ(5) = 2\big[\,v₂C(b−1,3) + R(5)\,\big].`

**The `j=5` tie.** `b−1` odd, so `C(b−1,3)` odd iff `b−1 ≡ 3 (mod 4)` iff `b ≡ 0 (mod 4)` (Lemma
B); thus `v₂C(b−1,3) = 0 ⟺ b ≡ 0 (mod 4)`. For `R(5) = v₂C(a+2,5) + v₂(K−20) − v₂(a+1)`: with
`A = v₂(a+1)`, `v₂C(a+2,5) = v₂C(a+1,4)`, and a direct computation of `v₂C(a+1,4)` (writing
`a+1 = 2^A w`, `w` odd) gives
`R(5) = 0 ⟺ A ≥ 2 ⟺ a ≡ 3 (mod 4)` (when `A=1`, `R(5) = v₂(\tfrac{a+1}{2}−1) ≥ 1 > 0`).
Therefore `\tilde Δ(5) = 0` iff **`a ≡ 3` and `b ≡ 0 (mod 4)`**, matching §0. ∎ (interior)

---

## 6. The two boundary-tie families

For `a` even the tie point `j=2` is interior iff `b ≥ 3`; for `a` odd the tie point `j=5` is
interior iff `b ≥ 6`. The two remaining tie carriers sit on the boundary `j = b+1` and are proved
directly from `M_{b+1} = b`, `M_b = 2b(m−b)` (Lemma 2):

**`b = 1`, `a` even** (`λ = (2m−2,1,1)`). `M₀ = (m−1)(2m−1)`, `M₁ = 2(m−1)`, `M₂ = 1`. Then
`val(0) = 2v₂(m−1)`, `val(2) = 2 + 2v₂C(m,2) = 2v₂(m) + 2v₂(m−1)`, so
`val(2) − val(0) = 2v₂(m) ≥ 0`, **`= 0` iff `m` odd**, i.e. `a = 2(m−1) ≡ 0 (mod 4)` — matching
`a ≡ 0, b = 1 ≡ 1 (mod 4)`. And `val(1) = 3 + 2v₂(m) + 2v₂(m−1) > val(0)`. ✓

**`b = 4`, `a` odd** (`λ = (2m−5,4,1)`). `M₃ = (a−3)(2a−1)`, `M₄ = 8(m−4)`, `M₅ = 4`. Then
`val(5) − val(3) = 2 v₂(m−3)` (after the falling-factorial reduction `C(m,5)/C(m,3) =
(m−3)(m−4)/20`), **`= 0` iff `m` even**, i.e. `a = 2m−5 ≡ 3 (mod 4)` — matching `a ≡ 3,
b = 4 ≡ 0 (mod 4)`; and `val(4) − val(3) = 1 + 2v₂(m−3) > 0`. ✓

*(Both families verified `m ≤ 200`.)*

---

## 7. The one residual gap (stated precisely) and verification

By parity (`val(j) ≡ j (mod 2)`) the **only** boundary point of the right parity to be an extra
minimizer is `j = b+1`: for `a` even `b+1` is even (matching `{0,2}`), for `a` odd `b+1` is odd
(matching `{3,5}`). For `b = 2` (`a` odd) `b+1 = 3 = j₀` is the genuine minimizer; the tie cases are
`b ∈ {1,4}` (§6). For **all other `b`** we need `val(b+1) > val(j₀)`, i.e. that `j=b+1` adds no
spurious minimizer (which would otherwise threaten the evenness). Unlike the interior `Δ`, this
comparison crosses between the factored form (`M_0/M_3`) and the boundary value `M_{b+1}=b`, and we
do **not** have a hand proof; it is verified **`m ≤ 80`** (`b ∉ {1,2,4}`), minimum margin
`val(b+1) − val(j₀) = 2`. This is the only step of the `c=1` theorem not established in general.

**Computational verification (all in `projects/code/threerow-c1/`):**
- Prop. 1 (general `M_j`) and Lemma 2 (`c=1`, all four forms + factorization of `B`) vs the
  symmetric-function/Murnaghan–Nakayama definition: all three-row `λ`, `m ≤ 12`. ✓
- Prop. 2 (`Δ(j)` Kummer formula) vs direct `val(j)−val(0)`: all `λ=(a,b,1)`, `m ≤ 15`. ✓
- Lemma C (compensation) and `v₂C(n,k) ≥ v₂(n)−v₂(k)`: all `λ=(a,b,1)`, `1 ≤ j ≤ b−1`, `m ≤ 30`. ✓
- Full Theorem (`J*` box + tie congruences): all `λ=(a,b,1)`, **`m ≤ 80`**, 0 violations
  (1170 ties, all `|J*| = 2`). ✓
- Boundary families `b=1, b=4` hand-formulas: `m ≤ 200`. ✓

---

## 8. Diagnosis — where the two-row collapse fails, and the general three-row picture

**The wall, named.** In two rows `M_j = f^{(a−j,b−j)}` is a *single* dimension and `val(j)−val(0)`
is a pure digit-sum expression with one bare `s₂(j)` (Lemma A closes it, parity `a≡b` killing
`j=1,3`). In three rows `M_j` is the signed `S₃`-determinant of Prop. 1. For `c=1` it collapses to
**three** binomials, and the Prop-2 reduction produces, beyond the two-row skeleton, the single new
pair `v₂(b(a+1) − j(j−1)) − v₂(a+1)` — **the valuation of a quadratic in `j`, not of a
factorial.** This is the concrete face of the `e₂ mod 2` wall flagged in both prior notes. What
*tames* it here is **Lemma C**: the quadratic-valuation term, added to `v₂C(a+2,j)`, always pays
back `v₂(a+1)` (`R(j) ≥ 0`). Lemma C is the three-row replacement for the two-row parity argument.

**The `K`-parity dichotomy.** Because `a+b = 2m−1` is odd, the parity of `a` decides everything:
`a` even ⟹ `K = b(a+1)` odd ⟹ the new term vanishes ⟹ the offset stays `j₀ = 0` and `J* ⊆ {0,2}`
(the two-row picture survives); `a` odd ⟹ `K ≡ 0 (mod 4)` ⟹ a **forced descent** `val(j) =
val(0)−j` for `j = 0,1,2,3` pushes the offset to `j₀ = 3` and `J* ⊆ {3,5}`. The tie at `j₀+2` is
governed in both cases by Lemma B (`C(·,2)`/`C(·,3)` odd) plus, for `a` odd, the exact value
`R(5)=0 ⟺ a ≡ 3 (mod 4)`.

**General three rows (`c` arbitrary).** Census (prior CODE + this session, `m ≤ 12`): three-row
`|J*| ∈ {1,2,4}`, always even when `≥ 2`. The new value `|J*| = 4` appears (e.g. `λ=(9,6,3)`,
`m=9`, `J* = {3,5,7,9} = 3 + 2·\{0,1,2,3\}`, an affine 2-adic box), so the strong `J* ⊆ {j₀,j₀+2}`
of the two-row/`c=1` families is **false for `c ≥ 2`**; the correct general statement is the
affine-2-adic-box conjecture (`|J*| = 2^{|S|}`). The `c=1` proof shows the mechanism that forces
evenness is *not* an argmin-internal symmetry but the arithmetic of `M_j`: specifically Lemma C
(compensation) controlling the quadratic-valuation term. For `c ≥ 2` the determinant of Prop. 1
keeps all six `S₃` terms and the analogue of Lemma 2 is a genuine signed multi-binomial; the
leading 2-adic layer is then a difference of several equal-valuation binomials, and "`|J*|` even"
becomes "this difference has an even number of surviving terms mod 2" — the precise open kernel.

---

### Files
`projects/code/threerow-c1/`: `mn.py` (Murnaghan–Nakayama), `dj.py` (Prop. 1), `c1d.py`/`prop2b.py`
(Lemma 2, Prop. 2), `general.py` (`Δ(j)`), `finalcheck.py` (Lemma C + box + ties, `m ≤ 80`),
`checkbdry.py` (boundary families).
