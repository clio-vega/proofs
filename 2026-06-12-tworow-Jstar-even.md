# Two-row shapes: `J* ⊆ {0,2}`, so `|J*|` is even on every tie (d = 4)

**Date:** 2026‑06‑12 (prove session)
**Status:** Two‑row family **PROVED** (complete), sharp form `J* ⊆ {0,2}`, ties characterised.

Companion to `2026-06-12-hook-Jstar-even.md`. Same setup, same conclusion shape
(`J* ⊆ {0,2}`), reached by the same skeleton: a **closed form for `M_j`**, then **honest
Kummer arithmetic** on `g(j)=val(j)−val(0)`. Two-row + hook together now cover two infinite
families of the general‑λ d=4 tie program.

---

## 0. The problem

Fix `m ≥ 1` and `λ ⊢ 2m`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`,

> `M_j = ⟨s_λ , h₁^{2(m−j)} e₂^j⟩ ∈ ℤ_{≥0}`,  `M₀ = f^λ`,
> `val(j) = j + 2 v₂( C(m,j) M_j )` (`= ∞` when `C(m,j)M_j = 0`),
> `V = min_j val(j)`,  `J* = { j : val(j) = V }`.

Leading‑π dichotomy (proved earlier): with `π = 1+i`,
`G_λ(i) = Σ_j C(m,j) M_j (i−1)^j` and `C ≡ |J*| (mod π)`, so **`|J*|` odd ⟹ `G_λ(i) ≠ 0`**;
the residual ties are exactly the shapes with `|J*|` even.

**This note proves**, for every two‑row `λ = (a,b) ⊢ 2m` (`a ≥ b ≥ 0`, `a+b = 2m`):

> **Theorem.** `0 ∈ J*` and `J* ⊆ {0,2}`. Hence `|J*| ∈ {1,2}`, so **`|J*|` is even whenever
> `|J*| ≥ 2`** (the leading‑π layer cancels on every two‑row tie). The unique possible tie is
> `J* = {0,2}`, occurring **iff `b ≡ 2,3 (mod 4)` and `a ≡ 1,2 (mod 4)`** (equivalently iff both
> `C(b,2)` and `C(a+1,2)` are odd).

Throughout, `s₂(n)` is the binary digit sum, `v₂` the 2‑adic valuation, `v₂(n!) = n − s₂(n)`
(Legendre), and `v₂ C(N,k) = s₂(k) + s₂(N−k) − s₂(N)` (Kummer).

---

## 1. The two‑row closed form

> **Lemma 1.** For `λ = (a,b) ⊢ 2m` (`a ≥ b ≥ 0`) and all `j ≥ 0`,
> `  M_j = C(2(m−j), b−j) − C(2(m−j), b−j−1) = f^{(a−j,\,b−j)},`
> the dimension of the two‑row shape with both rows shortened by `j` (so `M_j = 0` for `j > b`).
> Equivalently, the generating identity
> `  Σ_j C(m,j) M_j x^j = ([u^a] − [u^{a+1}])\,\big(u² + (2+x)u + 1\big)^m
>                      = ([u^b] − [u^{b−1}])\,\big(u² + (2+x)u + 1\big)^m.`

**Proof.** Since `h₁² = h₂ + e₂` (i.e. `s₂ + 2s_{11} = ...`; concretely `h₁²−e₂ = h₂`),
`h₁² + x e₂ = h₂ + (1+x)e₂`. Put `F = (h₁² + x e₂)^m`, a symmetric function of degree `2m`; then
`Σ_j C(m,j) M_j x^j = ⟨s_λ, F⟩`.

*Step 1 — reproducing kernel in two variables.* For the Hall inner product and any partition
`μ = (p,q)` with `p ≥ q`, `⟨F, h_p h_q⟩` equals the coefficient `c_{(p,q)}` of `m_{(p,q)}` in the
monomial expansion of `F`, which (restricting to two alphabet variables `s,t`) is the coefficient
`[s^p t^q] F(s,t,0,0,…)`. In two variables `h₁ = s+t` and `e₂ = st`, so
`  F(s,t) = \big((s+t)² + x\,st\big)^m = \big(s² + (2+x)\,st + t²\big)^m.`
`F` is homogeneous of degree `2m`; writing `u = s/t` and `g(u) := (u² + (2+x)u + 1)^m`,
`F(s,t) = t^{2m} g(u)`, hence for `p + q = 2m`, `[s^p t^q] F = [u^p] g =: g_p`.

*Step 2 — two‑row Jacobi–Trudi.* `s_{(a,b)} = h_a h_b − h_{a+1} h_{b−1}`, so (using `a ≥ b` and
`a+1 ≥ b−1`)
`  ⟨s_λ, F⟩ = [s^a t^b]F − [s^{a+1} t^{b−1}]F = g_a − g_{a+1}.`
Since `g` is palindromic (`u²+(2+x)u+1` is), `g_a = g_{2m−a} = g_b` and `g_{a+1} = g_{b−1}`, giving
the second form `g_b − g_{b−1}`.

*Step 3 — extract `[x^j]`.* Write `u²+(2+x)u+1 = (1+u)² + x u`, so by the binomial theorem
`  g_k = [u^k]\big((1+u)² + xu\big)^m = Σ_q C(m,q)\,x^q\,[u^{k−q}](1+u)^{2(m−q)}
       = Σ_q C(m,q)\,x^q\,C(2(m−q),\,k−q).`
Thus `[x^j] g_k = C(m,j)\,C(2(m−j),\,k−j)`. With `k = b` and `k = b−1`,
`  C(m,j) M_j = [x^j](g_b − g_{b−1}) = C(m,j)\big[C(2(m−j),b−j) − C(2(m−j),b−j−1)\big],`
and the ballot identity `C(N,r)−C(N,r−1) = f^{(N−r,\,r)}` with `N = 2(m−j)`, `r = b−j` gives the
`f^{(a−j,b−j)}` form (`N−r = a−j`). ∎

*Sanity.* `M₀ = C(2m,b) − C(2m,b−1) = f^{(a,b)}`, the standard two‑row dimension; setting `x = i−1`
(so `2+x = π = 1+i`) recovers the memory's two‑row engine
`G_λ(i) = [u^b](u²+πu+1)^m − [u^{b−1}](u²+πu+1)^m` (central‑trinomial / parameter‑Gegenbauer
object). Verified against the symmetric‑function definition via Murnaghan–Nakayama for all two‑row
`λ`, `m ≤ 8`. ✓

---

## 2. The Prop‑2 analogue: `g(j)` is a pure Kummer expression

For `0 ≤ j ≤ b` (where `M_j > 0`), the ballot/factorial form of Lemma 1 reads
`  C(m,j) M_j = (a−b+1)·\dfrac{m!\,(2m−2j)!}{j!\,(m−j)!\,(b−j)!\,(a−j+1)!}`
(from `M_j = C(2(m−j),b−j)·(a−b+1)/(a−j+1)`, the ballot form of the difference, and
`(a−j)!(a−j+1) = (a−j+1)!`).

> **Proposition 2 (two‑row).** For `0 ≤ j ≤ b`,
> `  g(j) := val(j) − val(0) = j + 2\big[\,v₂C(b,j) + v₂C(a+1,j) − s₂(j)\,\big].`

**Proof.** The constant factor `(a−b+1)` cancels in `g(j)`. Using `v₂(n!) = n − s₂(n)`,
`  v₂(C(m,j)M_j) − v₂(M_0) = N_{lin}(j) − S(j),`
where `N_{lin}` collects the linear (`n`) parts and `S` the digit‑sum parts. The linear parts cancel
identically:
`  N_{lin} = \big[m + (2m−2j) − j − (m−j) − (b−j) − (a−j+1)\big] − \big[2m − b − (a+1)\big] = 0`
(each of `m, a, b, j` and the constant has total coefficient `0`; uses `a+b = 2m`). The digit‑sum
part is
`  −S(j) = −\big[s₂(m) + s₂(2m−2j) − s₂(j) − s₂(m−j) − s₂(b−j) − s₂(a−j+1) − s₂(2m) + s₂(b) + s₂(a+1)\big].`
Now `s₂(2m) = s₂(m)` and `s₂(2m−2j) = s₂(m−j)` cancel in pairs, leaving
`  −S(j) = s₂(j) + s₂(b−j) + s₂(a+1−j) − s₂(b) − s₂(a+1).`
By Kummer, `s₂(j)+s₂(b−j)−s₂(b) = v₂C(b,j)` and `s₂(j)+s₂(a+1−j)−s₂(a+1) = v₂C(a+1,j)`, so
`−S(j) = v₂C(b,j) + v₂C(a+1,j) − s₂(j)`. Finally `g(j) = j + 2(N_{lin}−S(j)) = j + 2(−S(j))`. ∎

This is the exact two‑row counterpart of the hook Prop. 2: there the cancellation produced
`v₂C(m,j)+v₂C(b,j)−v₂C(2m−1,j)`; here it produces `v₂C(b,j)+v₂C(a+1,j)−s₂(j)`. The single bare
`s₂(j)` is the whole story.

---

## 3. The arithmetic lemmas

> **Lemma A (digit‑sum envelope).** For all `j ≥ 0`, `2 s₂(j) ≤ j + 1`, with equality **iff**
> `j ∈ {1,3}`. Moreover `2 s₂(j) ≤ j − 1` for every `j ≥ 4` (so `2s₂(j) − j ≤ −1` there), while
> `2s₂(j) = j` exactly for `j ∈ {0,2}`.

**Proof.** If `s₂(j) = k` then `j ≥ 2^k − 1` (the `k` ones sit in the lowest positions), so
`j + 1 ≥ 2^k`. Since `2^k ≥ 2k` for `k ≥ 1` (equality iff `k ∈ {1,2}`) and `2^0 = 1 > 0 = 2·0`,
`  2 s₂(j) = 2k ≤ 2^k ≤ j + 1,`
with equality only when `2k = 2^k` (`k ∈ {1,2}`) **and** `j = 2^k − 1` (`j ∈ {1,3}`). For `j ≥ 4`:
if `k ≥ 3` then `2^k ≥ 2k + 2`, so `j+1 ≥ 2^k ≥ 2k+2 = 2s₂(j)+2`; if `k ≤ 2` and `j ≥ 4` then
`j ≥ 5` when `k = 2` and `j ≥ 4` when `k = 1`, so `2s₂(j) = 2k ≤ 4 ≤ j−1`. Either way
`2s₂(j) ≤ j−1`. The case `2s₂(j) = j` forces `j` even and (writing `j = 2i`) `s₂(i) = i`, i.e.
`i ∈ {0,1}`, so `j ∈ {0,2}`. ∎

> **Lemma B (small binomials mod 2).** `C(n,2)` is odd `⟺ n ≡ 2,3 (mod 4)`; `C(n,3)` is odd
> `⟺ n ≡ 3 (mod 4)`.

**Proof.** `v₂C(n,2) = v₂(n(n−1)) − 1`; of the two consecutive integers exactly one is even, and
`C(n,2)` is odd iff that even one is `≡ 2 (mod 4)`, i.e. `n ≡ 2,3 (mod 4)`. For `C(n,3) =
n(n−1)(n−2)/6`, `v₂C(n,3) = v₂(n(n−1)(n−2)) − 1`. If `n` is odd, only `n−1` is even and
`C(n,3)` is odd iff `v₂(n−1) = 1`, i.e. `n ≡ 3 (mod 4)`. If `n` is even, both `n` and `n−2` are
even and exactly one is `≡ 0 (mod 4)`, so `v₂(n(n−2)) ≥ 3` and `C(n,3)` is even. ∎

---

## 4. Proof of the Theorem

Write `W(j) := v₂C(b,j) + v₂C(a+1,j) ≥ 0`, so `g(j) = j + 2W(j) − 2s₂(j)` (Prop. 2). Note
`g(0) = 0`. We show `g(j) ≥ 0` for all `0 ≤ j ≤ b`, with `g(j) > 0` for `j ∉ {0,2}`; then `val(0)`
is the minimum (`0 ∈ J*`) and `J* ⊆ {0,2}`. Crucially, `a + b = 2m` is even, so **`a ≡ b (mod 2)`**.

**`j ≥ 4`.** By Lemma A, `2s₂(j) ≤ j − 1`, so `g(j) ≥ j + 0 − (j−1) = 1 > 0`. (Unconditional.)

**`j = 1`** (needs `b ≥ 1`). `2s₂(1) = 2`, so `g(1) = 2W(1) − 1` with `W(1) = v₂(b) + v₂(a+1)`. If
`b` is even, `v₂(b) ≥ 1`; if `b` is odd then `a` is odd, so `a+1` is even and `v₂(a+1) ≥ 1`. Either
way `W(1) ≥ 1`, hence `g(1) ≥ 1 > 0`.

**`j = 3`** (needs `b ≥ 3`). `2s₂(3) = 4`, so `g(3) = 2W(3) − 1` with
`W(3) = v₂C(b,3) + v₂C(a+1,3)`. By Lemma B, `C(b,3)` is odd `⟺ b ≡ 3 (mod 4)` (forcing `b` odd),
and `C(a+1,3)` is odd `⟺ a+1 ≡ 3 (mod 4) ⟺ a ≡ 2 (mod 4)` (forcing `a` even). Both odd would need
`b` odd and `a` even, contradicting `a ≡ b (mod 2)`. So at least one is even, `W(3) ≥ 1`, and
`g(3) ≥ 1 > 0`.

**`j = 2`** (needs `b ≥ 2`). `2s₂(2) = 2`, so `g(2) = 2W(2) = 2[v₂C(b,2) + v₂C(a+1,2)] ≥ 0`, with
`g(2) = 0 ⟺ C(b,2)` and `C(a+1,2)` both odd. By Lemma B this is `b ≡ 2,3 (mod 4)` and
`a+1 ≡ 2,3 (mod 4)`, i.e. `a ≡ 1,2 (mod 4)`.

Thus `g(j) > 0` for every `j ∉ {0,2}` and `g(2) ≥ 0`, so `J* ⊆ {0,2}`, `0 ∈ J*`, and the tie
`J* = {0,2}` occurs exactly when `b ≡ 2,3 (mod 4)` and `a ≡ 1,2 (mod 4)`. In every tie
`|J*| = 2` is even. ∎

**Remark (vanishing).** Combined with the dichotomy, a two‑row `G_λ(i) = 0` forces `|J*|` even,
hence `J* = {0,2}`, hence `b ≡ 2,3 (mod 4)` and `a ≡ 1,2 (mod 4)` — pinning the vanishing
candidates to that congruence box and matching the `b ≡ 2,3 (mod 4)` frontier recorded in memory.
The unique known two‑row vanisher `(2,2)` (`a=b=2`, `m=2`) satisfies it: `b = 2 ≡ 2`, `a = 2 ≡ 2`. ✓

---

## 5. Verification (computational)

- **Lemma 1**, `M_j = C(2(m−j),b−j) − C(2(m−j),b−j−1)`, against the symmetric‑function definition
  `⟨s_λ,h₁^{2(m−j)}e₂^j⟩` computed by Murnaghan–Nakayama (`M_j = 2^{−j}Σ_t C(j,t)(−1)^t
  χ^λ(2^t 1^{2(m−t)})`): all two‑row `λ`, `m ≤ 8`. ✓
  The boxed GF identity `Σ C(m,j)M_j x^j = ([u^a]−[u^{a+1}])(u²+(2+x)u+1)^m`: symbolic in `x`,
  all two‑row `λ`, `m ≤ 7`. ✓
- **Proposition 2** `g(j) = j + 2[v₂C(b,j)+v₂C(a+1,j)−s₂(j)]` against direct `val(j)−val(0)`:
  all two‑row `λ`, `m ≤ 14`, 0 violations. ✓
- **Theorem**: `J* ⊆ {0,2}`, `0 ∈ J*`, and tie `⟺ C(b,2),C(a+1,2)` both odd `⟺ b ≡ 2,3 (mod 4)
  and a ≡ 1,2 (mod 4)`: all two‑row `λ`, `m ≤ 14`, 0 violations (28 ties found). ✓
- **Lemmas A, B**: `2s₂(j) ≤ j+1` (eq. iff `j ∈ {1,3}`), `2s₂(j) ≤ j−1` for `j ≥ 4`;
  `C(n,2)` odd `⟺ n ≡ 2,3 (mod 4)`, `C(n,3)` odd `⟺ n ≡ 3 (mod 4)`: checked `n,j < 2·10^5`. ✓

---

## 6. Where this leaves the general program

Two infinite families now closed by the **same template**:

| family | closed form `M_j` | Prop‑2 `g(j)` | result |
|---|---|---|---|
| hook `(a,1^b)` | `C(2m−1−j, a−1)` | `j − 2(c₃−c₁−c₂)` | `J* ⊆ {0,2}` |
| two‑row `(a,b)` | `C(2(m−j),b−j) − C(2(m−j),b−j−1) = f^{(a−j,b−j)}` | `j + 2[v₂C(b,j)+v₂C(a+1,j)−s₂(j)]` | `J* ⊆ {0,2}` |

Both reduce the tie to a **single bare digit‑sum envelope** (`s₂` of `t` for hooks, of `j` for
two‑rows) controlled by Lemma A, with the parity constraint `a ≡ b (mod 2)` killing the only escape
points `j = 1, 3`. In both families the answer is the *strong* form `J* ⊆ {0,2}`, not merely
`|J*|` even.

**The next obstruction.** For three or more rows, the reproducing‑kernel + Jacobi–Trudi step still
gives `Σ C(m,j)M_j x^j = ([s^{α}t^{β}…] \text{ of a power})`, but `M_j` is then a *signed* sum of
several binomials (the Jacobi–Trudi determinant has `> 2` terms), so the falling‑factorial
collapse of Prop. 2 no longer telescopes to one digit‑sum term. The precise general statement
(§6.2 of the hook note: the slope‑`(−1/2)` residual polynomial `e(y)` is a power of `(1+y)`)
remains open, but the two‑row case confirms its mechanism: it is the parity constraint
`a ≡ b (mod 2)` — a *symmetry of the shape*, not of the Newton polygon — that forces evenness, as
predicted in §1.3 of the hook note.

---

### Files
`code:` `/tmp/tworow1.py` (GF↔closed‑form), `/tmp/tworow_mn.py` (M_j via Murnaghan–Nakayama),
`/tmp/tworow_g.py` (g(j), J*, tie char), `/tmp/lemmas.py` (Lemmas A,B). To migrate into
`projects/code/`.
