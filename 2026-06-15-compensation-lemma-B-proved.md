# Compensation Lemma B (the `a`-odd two-generator Number Lemma), three-row `c=3` — PROVED

**Date:** 2026-06-15 (prove session)
**Status:** **FULL WIN.** The two-generator inequality `Δ̃(j) ≥ 0` (Gap 2, the `e₂ mod 2` wall in its
first concrete two-generator instance) is proved unconditionally for all interior `4 ≤ j ≤ b`, together
with the complete tie classification: equality holds **only** at `j ∈ {5,7,9}`, all three tie together,
and they tie **iff `a ≡ 1, b ≡ 2 (mod 4)`**. Hence in the `a`-odd branch `J*` is `{3}` or the full box
`{3,5,7,9} = 3 + ⟨2,4⟩`, so `|J*| ∈ {1,4}` — closing the `c=3` interior modulo the standard boundary
residual (Gap 3).

Companion to and completion of `2026-06-14-threerow-c3-Jstar-even.md` (§5, §7 Gap 2). Builds on the
proved general-`c` Number Lemma (`2026-06-14-numberlemma-general-c.md`); in fact the whole proof bottoms
out at the elementary digit-sum bound `2 s₂(k) ≤ k+1`, exactly the "Gap-1 doubled" shape predicted.

---

## 0. Statement

Fix `λ = (a,b,3) ⊢ 2m` with `a` **odd**, `b` **even**, `a ≥ b ≥ 4`, `a+b = 2m−3`. After the proved
forced descent (Prop. 3 of the `c=3` note), the offset is `j₀ = 3` and `Δ̃(j) := val(j) − val(3)`.

> **Compensation Lemma B.** For all interior `4 ≤ j ≤ b`,
> `Δ̃(j) = j + 3 − 2 s₂(j) + 2 U(j) ≥ 0`,  where
> `U(j) := v₂C(a+4,j) + v₂C(b+3,j) + v₂Q₃(j) − v₂Q₃(0)`,
> with **equality only at `j ∈ {5,7,9}`**, in which case
> `Δ̃(5)=Δ̃(7)=Δ̃(9)=0 ⟺ a ≡ 1 and b ≡ 2 (mod 4)`.

Notation: `s₂` = binary digit sum, `v₂` = 2-adic valuation, `v₂C(n,k) = s₂(k)+s₂(n−k)−s₂(n)` (Kummer),
`C(n,k)=0` for `k<0` or `k>n`. Falling factorial `n^{\underline{k}} = n(n−1)\cdots(n−k+1) = k!\,C(n,k)`.

**Inputs used (proved earlier, cited, not re-derived):**
- *Closed form / Prop. 2:* `U(j) = v₂P(j) − v₂P(0)` with `P(j) := C(a+4,j)\,C(b+3,j)\,Q₃(j)`
  (`P(0)=Q₃(0)`), and `Δ̃(j) = j+3−2s₂(j)+2U(j)`.
- *Structural identity:* `Q₃(j) = (a−1)(b−2)\,H(j) − 720\,C(j,6)`, with `720 = 6!`, and
  `H(j) = (a+3)(b+2)\,G(j) − 6\,E(j)`, where
  `G(j) = (a+4)(b+3) − 6j`,  `E(j) = C(j,2)\,Φ(j)`,  `Φ(j) = ab+a+2b − (j−1)(j−2)`.
  *(The `H` split and `E = C(j,2)Φ` are algebraic identities, machine-checked symbolically; see Files.)*
- *Forced descent (Prop. 3):* `j₀ = 3`, so `Δ̃ = Δ + 3`.

**Standing notation for valuations** (all `≥ 1`, since `a` odd, `b` even):
`α = v₂(a−1)`, `α' = v₂(a+3)`, `β = v₂(b−2)`, `β' = v₂(b+2)`.
Then `a, a+2, a+4` and `b±1, b±3` are odd, and
`v₂Q₃(0) = v₂[(a−1)(a+3)(a+4)(b−2)(b+2)(b+3)] = α+α'+β+β'`,
`v₂H(0)  = v₂[(a+3)(a+4)(b+2)(b+3)] = α'+β'`.

---

## 1. Three elementary lemmas

Set `g(j) := s₂(j) − (j+3)/2`. Since `val(j) ≡ j (mod 2)`, `Δ̃(j) = 2(U(j) − g(j))`, so

> **`Δ̃(j) ≥ 0 ⟺ U(j) ≥ g(j) ⟺ v₂P(j) ≥ v₂P(0) + g(j)`.**   (the *target* T(j) := v₂P(0)+g(j))

**Lemma A.** `2 s₂(k) ≤ k+1` for all `k ≥ 0`, with equality iff `k ∈ {1,3}`.

*Proof.* Induction. `k=0`: `0≤1`. For `k ≥ 1` write `k = 2q+r`, `r ∈ {0,1}`; then `s₂(k)=s₂(q)+r`.
By IH `2s₂(q) ≤ q+1`, so `2s₂(k) = 2s₂(q)+2r ≤ q+1+2r`, and `q+1+2r ≤ 2q+r+1 = k+1 ⟺ r ≤ q`. If
`r=0` this is `0 ≤ q`; if `r=1` it is `q ≥ 1`, i.e. `k ≥ 3`, and the lone remaining case `k=1` gives
`2 ≤ 2` directly. Equality forces `r=q` at every step: `k=1` (`q=0,r=1`) and `k=3` (`q=1,r=1`). ∎

**Lemma C.** `2 s₂(j) + 2 v₂(j(j−1)) ≤ j+3` for all `j ≥ 2`.

*Proof.* Exactly one of `j, j−1` is even.
*Case `j` even,* `j = 2^v q` (`q` odd, `v ≥ 1`), `v₂(j−1)=0`, `s₂(j)=s₂(q)`. Need
`2s₂(q)+2v ≤ 2^v q+3`. By Lemma A `2s₂(q) ≤ q+1`, so it suffices that `2v ≤ q(2^v−1)+2`; since
`q ≥ 1` and `2^v+1 ≥ 2v` for `v ≥ 1`, `q(2^v−1)+2 ≥ 2^v+1 ≥ 2v`.
*Case `j` odd,* `j−1 = 2^w r` (`r` odd, `w ≥ 1`), `v₂(j)=0`, `s₂(j)=s₂(r)+1`. Need
`2s₂(r)+2+2w ≤ 2^w r+4`, i.e. `2w ≤ r(2^w−1)+1`; since `r ≥ 1`, `r(2^w−1)+1 ≥ 2^w ≥ 2w`. ∎

**Lemma NL₁** (the `c=1` Number Lemma, self-contained). For `F` even `≥ 2` and `2 ≤ j ≤ F+1`,
`v₂C(F+1,j) + v₂(j(j−1)) ≥ v₂(F)`.

*Proof.* The subset identity `C(F+1,j)\,C(j,2) = C(F+1,2)\,C(F−1,j−2)` gives
`v₂C(F+1,j) + v₂C(j,2) = v₂C(F+1,2) + v₂C(F−1,j−2)`. Now `2\,C(F+1,2) = (F+1)F`, so
`v₂C(F+1,2) = v₂(F) − 1` (`F+1` odd); and `2\,C(j,2) = j(j−1)`, so `v₂C(j,2) = v₂(j(j−1)) − 1`.
Substituting and using `v₂C(F−1,j−2) ≥ 0`,
`v₂C(F+1,j) + v₂(j(j−1)) = v₂(F) + v₂C(F−1,j−2) ≥ v₂(F)`. ∎

Two consequences, used below.

**Lemma (i).** For `b` even and `4 ≤ j ≤ b`,  `v₂C(b+3,j) ≥ β' + g(j)`.
*Proof.* NL₁ with `F = b+2` (even): `v₂C(b+3,j) ≥ v₂(b+2) − v₂(j(j−1)) = β' − v₂(j(j−1))`. By Lemma C,
`v₂(j(j−1)) ≤ (j+3)/2 − s₂(j) = −g(j)`, so `v₂C(b+3,j) ≥ β' + g(j)`. ∎

**Lemma B.** For `a` odd and `6 ≤ j ≤ a`,  `v₂C(a+4,j) ≥ α + α' + g(j) − 1`.
*Proof.* The subset identity `C(a+4,j)\,C(j,6) = C(a+4,6)\,C(a−2,j−6)` gives
`v₂C(a+4,j) + v₂C(j,6) = v₂C(a+4,6) + v₂C(a−2,j−6)`. Here `720\,C(a+4,6) = (a+4)^{\underline 6}` whose
three even factors are `a+3,a+1,a−1`, so `v₂C(a+4,6) = (α+α'+v₂(a+1)) − 4`. Hence
`v₂C(a+4,j) = α+α'+v₂(a+1) − 4 − v₂C(j,6) + v₂C(a−2,j−6)`.
Plugging into the claim and cancelling `α+α'`, Lemma B is equivalent to
> **(B′)** `v₂(a+1) + v₂C(a−2,j−6) ≥ R(j)`,  `R(j) := v₂C(j,6) + s₂(j) − (j+3)/2 + 3`.

Now `R(j) ≤ 1`: by Legendre `v₂C(j,6) = v₂(j!) − v₂(6!) − v₂((j−6)!) = 2 + s₂(j−6) − s₂(j)`, so the
`s₂(j)` cancels and `R(j) = 5 + s₂(j−6) − (j+3)/2`, i.e. (clearing the half)
`2R(j) = 2 s₂(j−6) − (j−6) + 1 ≤ 2` by Lemma A applied to `k = j−6` (`2s₂(k) ≤ k+1`), with equality iff
`j−6 ∈ {1,3}`, i.e. `j ∈ {7,9}`. Hence `R(j) ≤ 1`. Since `v₂(a+1) ≥ 1` (`a` odd) and
`v₂C(a−2,j−6) ≥ 0`, (B′) holds. ∎

---

## 2. The inequality `Δ̃(j) ≥ 0`

Write `P(j) = P₁(j) − P₂(j)` with
`P₁ = (a−1)(b−2)\,C(a+4,j)\,C(b+3,j)\,H(j)` (the **heavy** part) and
`P₂ = 720\,C(a+4,j)\,C(b+3,j)\,C(j,6)` (the **tip**), from `Q₃ = (a−1)(b−2)H − 720C(j,6)`.

By the ultrametric inequality `v₂P ≥ min(v₂P₁, v₂P₂)`, so it suffices to show **each** of `v₂P₁`,
`v₂P₂` is `≥ T(j) = α+α'+β+β' + g(j)`.

### 2a. Tip: `v₂P₂ ≥ T`

If `j < 6` then `C(j,6)=0`, `P₂=0`, trivially true. Let `6 ≤ j ≤ b`. Absorb `C(j,6)` into the
`(b+3)`-binomial: `C(b+3,j)\,C(j,6) = C(b+3,6)\,C(b−3,j−6)`, and `720\,C(b+3,6) = (b+3)^{\underline 6}`
whose three even factors are `b+2, b, b−2`, so `v₂ = β+β'+v₂(b)`. Therefore
`v₂P₂ = β+β'+v₂(b) + v₂C(a+4,j) + v₂C(b−3,j−6)`.
By Lemma B, `v₂C(a+4,j) ≥ α+α'+g(j)−1`, and `v₂(b) ≥ 1`, `v₂C(b−3,j−6) ≥ 0`; hence
`v₂P₂ ≥ β+β' + 1 + (α+α'+g(j)−1) = α+α'+β+β'+g(j) = T(j)`. ✓

### 2b. Heavy: `v₂P₁ ≥ T`

`v₂P₁ = (α+β) + v₂[\,C(a+4,j)\,C(b+3,j)\,H(j)\,]`, so it suffices that
`v₂[C(a+4,j)C(b+3,j)H(j)] ≥ α'+β'+g(j)`. Split `H = (a+3)(b+2)G − 6E`, giving
`C(a+4,j)C(b+3,j)H = Π₁ − Π₂`:

- **`Π₁ = (a+3)(b+2)·\big[C(a+4,j)C(b+3,j)G(j)\big]`.** The bracket is an integer, so
  `v₂Π₁ ≥ v₂((a+3)(b+2)) = α'+β' ≥ α'+β'+g(j)` because `g(j) ≤ 0`. ✓
- **`Π₂ = 6\,C(a+4,j)C(b+3,j)\,C(j,2)\,Φ(j)`.** Here `Φ = ab+a+2b − (j−1)(j−2)` is **odd**
  (`ab+a+2b` is odd since `a` odd, `b` even; `(j−1)(j−2)` is even), so `v₂Φ = 0`, and `v₂(6)=1`.
  Absorb `C(j,2)` into the `(a+4)`-binomial: `C(a+4,j)C(j,2) = C(a+4,2)C(a+2,j−2)`, and
  `2C(a+4,2) = (a+4)(a+3)` gives `v₂C(a+4,2) = α'−1`. Thus
  `v₂Π₂ = 1 + (α'−1) + v₂C(a+2,j−2) + v₂C(b+3,j) = α' + v₂C(a+2,j−2) + v₂C(b+3,j)`.
  By Lemma (i), `v₂C(b+3,j) ≥ β'+g(j)`, and `v₂C(a+2,j−2) ≥ 0`, so `v₂Π₂ ≥ α'+β'+g(j)`. ✓

Both `v₂Π₁, v₂Π₂ ≥ α'+β'+g(j)`, hence `v₂[Π₁−Π₂] ≥ α'+β'+g(j)`, i.e. `v₂P₁ ≥ T(j)`. ✓

Therefore `v₂P(j) ≥ T(j) = v₂P(0)+g(j)`, i.e. `U(j) ≥ g(j)`, i.e.
**`Δ̃(j) = 2(U(j)−g(j)) ≥ 0`** for all `4 ≤ j ≤ b`. ∎

---

## 3. The tie set is exactly `{5,7,9}`

`Δ̃(j) = 2(U(j)−g(j))`, so a tie at `j` means `U(j)=g(j)`, equivalently both `v₂P₁` and `v₂P₂` meet
their bound `T(j)` with no lift. We rule out all `j ∉ {5,7,9}`.

**Even `j`.** `Δ̃(j) = val(j) − val(3) ≡ j − 3 ≡ 1 (mod 2)`: odd. Being `≥ 0` it is `≥ 1`. No tie.

**Odd `j ≥ 11`: `Δ̃(j) ≥ 2`.** We show *both* pieces clear `T` with **slack `≥ 1`**, so `v₂P ≥ T+1`,
`U ≥ g+1`, `Δ̃ ≥ 2`.

- *Tip slack.* From 2a, `v₂P₂ − T = (v₂(b)−1) + v₂C(b−3,j−6) + [v₂C(a+4,j) − (α+α'+g−1)]`. The last
  bracket is the slack of Lemma B, which by the proof of (B′) equals
  `v₂(a+1) + v₂C(a−2,j−6) − R(j)`. For odd `j ≥ 11`, `R(j) ≤ 0` (since `2R(j) ≤ 2` with equality only
  at `j ∈ {7,9}`, Lemma B remark), so the bracket is `≥ v₂(a+1) ≥ 1`. Hence `v₂P₂ − T ≥ 1`.
- *Heavy slack.* `Π₁` has slack `≥ −g(j) ≥ 2` (Lemma A ⟹ `g(j) ≤ −2` for `j ≥ 4`). For `Π₂`, the
  slack is `≥` slack of Lemma (i) `= v₂C(b+3,j) − (β'+g) ≥ −v₂(j(j−1)) − g = −v₂(j−1) − g` (`j` odd).
  For odd `j ≥ 11`, `2s₂(j) + 2v₂(j−1) ≤ j+1` (a strict form of Lemma C — same case proof, with the
  small exceptional values `j=3,5,9` falling outside `j ≥ 11`), which rearranges to
  `−v₂(j−1) − g(j) ≥ 1`. So `min(v₂Π₁,v₂Π₂)` clears its bound by `≥ 1`, giving `v₂P₁ − T ≥ 1`.

So `Δ̃(j) ≥ 2` for odd `j ≥ 11`: no ties beyond the box.

**`j ∈ {5,7,9}`** are handled exactly in §4. Together with the above, the tie set is `⊆ {5,7,9}`.

---

## 4. Exact values at `j = 5, 7, 9` and the `(mod 4)` switch

The valuations of the relevant binomials are, by counting even factors in a window of consecutive
integers (a odd, b even):
`v₂C(a+4,5) = α'+v₂(a+1)−3`, `v₂C(b+3,5) = β'+v₂(b)−3`;
`v₂C(a+4,7) = α+α'+v₂(a+1)−4`, `v₂C(b+3,7) = β+β'+v₂(b)−4`;
`v₂C(a+4,9) = α+α'+v₂(a+1)+v₂(a−3)−7`, `v₂C(b+3,9) = β+β'+v₂(b)+v₂(b−4)−7`.
Recall the residue dictionary (a odd, b even):
`a≡1 (4) ⟺ α≥2, v₂(a+1)=1, v₂(a−3)=1`;  `a≡3 (4) ⟺ α=1, v₂(a+1)≥2, v₂(a−3)≥2`;
`b≡2 (4) ⟺ β≥2, v₂(b)=1, v₂(b−4)=1`;     `b≡0 (4) ⟺ β=1, v₂(b)≥2, v₂(b−4)≥2`.
Note also `α' ≥ 2 ⟺ a≡1`, `α'=1 ⟺ a≡3`; `β' ≥ 2 ⟺ b≡2`, `β'=1 ⟺ b≡0`.

### `j = 5` (tip vanishes: `C(5,6)=0`, so `Q₃(5)=(a−1)(b−2)H(5)`).
`Δ̃(5) = 8 − 4 + 2U(5) = 4 + 2U(5)`; tie `⟺ U(5) = −2`. Substituting the binomials and
`v₂Q₃(5)−v₂Q₃(0) = v₂H(5)−v₂H(0)`,
`U(5) = v₂(a+1) + v₂(b) + v₂H(5) − 6`.
From `H(5) = (a+3)(b+2)\,[(a+4)(b+3)−30] − 60\,(ab+a+2b−12)`: the first term has `v₂ = α'+β'` (odd
bracket), the second has `v₂ = 2` (`60 = 4·15`, `ab+a+2b−12` odd). Thus `v₂H(5)=2` whenever
`α'+β' > 2`; when `α'+β'=2` the two terms can only raise `v₂H(5) ≥ 2`.
- `(a,b)≡(1,2)`: `α'+β' ≥ 4 ⟹ v₂H(5)=2`; `v₂(a+1)=v₂(b)=1` ⟹ `U(5)=1+1+2−6=−2`: **tie**.
- `(1,0)`: `v₂H(5)=2`, `v₂(a+1)=1`, `v₂(b)≥2` ⟹ `U(5) ≥ −1`, `Δ̃(5) ≥ 2`.
- `(3,2)`: `v₂H(5)=2`, `v₂(a+1)≥2`, `v₂(b)=1` ⟹ `U(5) ≥ −1`, `Δ̃(5) ≥ 2`.
- `(3,0)`: `v₂(a+1)≥2`, `v₂(b)≥2`, `v₂H(5)≥2` ⟹ `U(5) ≥ 0`, `Δ̃(5) ≥ 4`.

### `j = 7`.
`Δ̃(7) = 10 − 6 + 2U(7) = 4 + 2U(7)`; tie `⟺ U(7)=−2`.
`U(7) = v₂(a+1) + v₂(b) + v₂Q₃(7) − 8`.
`H(7) = (a+3)(b+2)\,[(a+4)(b+3)−42] − 126\,(ab+a+2b−30)`: first term `v₂=α'+β' ≥ 2`, second
`v₂ = 1` (`126 = 2·63`, odd cofactor), so `v₂H(7) = 1` **always**. Then
`Q₃(7) = (a−1)(b−2)H(7) − 5040`, with `v₂` of the heavy part `= α+β+1` and `v₂(5040)=4`:
`v₂Q₃(7) = 3` if `α+β=2`; `= 4` if `α+β ≥ 4`; `≥ 4` if `α+β=3`.
- `(1,2)`: `α,β ≥ 2 ⟹ α+β ≥ 4 ⟹ v₂Q₃(7)=4`; `v₂(a+1)=v₂(b)=1` ⟹ `U(7)=1+1+4−8=−2`: **tie**.
- `(1,0)`: `α+β ≥ 3 ⟹ v₂Q₃(7) ≥ 4`; `v₂(a+1)=1`, `v₂(b)≥2` ⟹ `U(7) ≥ −1`, `Δ̃(7) ≥ 2`.
- `(3,2)`: `α+β ≥ 3 ⟹ v₂Q₃(7) ≥ 4`; `v₂(a+1)≥2`, `v₂(b)=1` ⟹ `U(7) ≥ −1`, `Δ̃(7) ≥ 2`.
- `(3,0)`: `α+β=2 ⟹ v₂Q₃(7)=3`; `v₂(a+1),v₂(b) ≥ 2` ⟹ `U(7) ≥ −1`, `Δ̃(7) ≥ 2`.

### `j = 9`.
`Δ̃(9) = 12 − 4 + 2U(9) = 8 + 2U(9)`; tie `⟺ U(9)=−4`.
`U(9) = v₂(a+1)+v₂(a−3)+v₂(b)+v₂(b−4) + v₂Q₃(9) − 14`.
`H(9) = (a+3)(b+2)\,[(a+4)(b+3)−54] − 216\,(ab+a+2b−56)`: first `v₂=α'+β'`, second `v₂=3`
(`216 = 8·27`), so `v₂H(9)=2` if `α'+β'=2`, `=3` if `α'+β' ≥ 4`, `≥3` if `α'+β'=3`. With
`Q₃(9) = (a−1)(b−2)H(9) − 60480`, `v₂(60480)=6`, heavy `v₂ = α+β+v₂H(9)`:
- `(1,2)`: `v₂H(9)=3`, heavy `v₂ = α+β+3 ≥ 7 > 6 ⟹ v₂Q₃(9)=6`; the four binomial valuations are
  `1,1,1,1` ⟹ `U(9) = 1+1+1+1+6−14 = −4`: **tie**.
- `(1,0)`: `v₂H(9) ≥ 3`, heavy `v₂ ≥ 6`, so `v₂Q₃(9) ≥ 6`; `v₂(b),v₂(b−4) ≥ 2` ⟹
  sum `≥ 1+1+2+2+6 = 12`, `U(9) ≥ −2`, `Δ̃(9) ≥ 4`.
- `(3,2)`: symmetric, sum `≥ 2+2+1+1+6 = 12`, `Δ̃(9) ≥ 4`.
- `(3,0)`: `α'+β'=2 ⟹ v₂H(9)=2`, heavy `v₂ = α+β+2 = 4 < 6 ⟹ v₂Q₃(9)=4`; the four binomial
  valuations are all `≥ 2` ⟹ sum `≥ 2+2+2+2+4 = 12`, `Δ̃(9) ≥ 4`.

**Conclusion.** For each `j ∈ {5,7,9}`, `Δ̃(j)=0 ⟺ (a,b) ≡ (1,2) (mod 4)`, and `Δ̃(j) ≥ 2` otherwise.
Combined with §2 (`Δ̃ ≥ 0`) and §3 (no ties for even `j` or odd `j ≥ 11`):

> `J* − 3` is `∅` if `(a,b) ≢ (1,2) (mod 4)`, and the full box `{2,4,6} = ⟨2,4⟩∖{0}`'s shift
> `{5,7,9}−3` if `(a,b) ≡ (1,2) (mod 4)`. Hence `J* = {3}` or `{3,5,7,9}`, so `|J*| ∈ {1,4}`. ∎

---

## 5. Why it works (the structure)

The whole proof is "**Gap-1 doubled**". The single `a`-even binomial split is replaced by a **two-sided**
one: the `(2c)!\,C(j,2c)` *tip* feeds its `2c` consecutive integers into **one** binomial to manufacture
that side's even-neighbour valuation, while the *heavy* factor's matching valuation is supplied by the
**other** binomial via the dual absorption — and the spare `α+α'` (resp. `β'`) **cancels**, collapsing
everything onto the elementary `2 s₂(k) ≤ k+1`. Concretely:

- The **tip** (generator `4`, tight at `j=7`) needs `α+α'+β+β'`; the `b`-side absorption of `C(j,6)`
  delivers `β+β'+v₂(b)` from `(b+3)^{\underline 6}`, and Lemma B delivers `α+α'` from the `a`-side, the
  residue cancelling down to `R(j) ≤ 1 ⟺ 2s₂(j−6) ≤ (j−6)+1`.
- The **heavy** (generator `2`, tight at `j=5`) factors as `(a+3)(b+2)G − 6\,C(j,2)Φ` with `Φ` odd;
  `Π₁` carries `α'+β'` outright, and `Π₂`'s `C(j,2)` absorbs into the `a`-side while Lemma (i)
  (`= NL₁ + Lemma C`) supplies `β'` from the `b`-side.

The two generators **lock** because both tight conditions (`v₂H(5)`, `v₂Q₃(7)`, `v₂Q₃(9)` minimal *and*
the binomial valuations minimal) are governed by the **same** pair of residues `a≡1, b≡2 (mod 4)`: each
makes `α',β' ≥ 2` (so the heavy/`H`-terms reach their minimal valuations) **and** `v₂(a+1)=v₂(b)=1` (so
the binomials reach theirs) **simultaneously**. That single coincidence is the `|J*|=4` event.

This is the first hand-proved **two-generator** Number Lemma. The engine is visibly `c`-graded
(`Q_{(a,b,c)} = (\text{heavy})·H − (2c)!\,C(j,2c)`, `H` carrying generators `2,…,2c−2`, the tip carrying
`2c`); the present proof is the template for the general `e₂ mod 2` wall.

---

## 6. Verification

Every link was machine-checked with **0 failures** over all `a`-odd shapes `(a,b,3)`, `4 ≤ j ≤ b`,
`m ≤ 79` (and the standalone digit lemmas to `k,j < 3·10⁴–2·10⁵`):
Lemma A, Lemma C, NL₁, Lemmas (i) & B, the heart `2R(j) ≤ 2`, the tip bound `v₂P₂ ≥ T`, the heavy bound
`v₂P₁ ≥ T`, the main inequality `Δ̃ ≥ 0`; the §4 binomial- and `U`-formulas and every per-residue case
split; the strictness `Δ̃(j) ≥ 2` for odd `j ≥ 11` and the two strictness ingredients `R(j) ≤ 0`,
`2s₂(j)+2v₂(j−1) ≤ j+1`. Symbolic identities `H = (a+3)(b+2)G − 6E` and `E = C(j,2)Φ` confirmed by
`sympy`.

## 7. What remains

Only **Gap 3**, the *standard boundary residual* `val(b+i) > val(j₀)` for the non-tie boundary points
`j ∈ {b+1,b+2,b+3}` (where the factored interior form gives way to boundary values) — verified `m ≤ 79`,
identical in nature to the lone residual in the `c=1,2` families. The entire `a`-odd **interior** is now
proved unconditionally.

### Files
`projects/code/threerow-c3/`: `c3gap2_explore.py`, `c3gap2_split.py`, `c3gap2_tip.py`,
`c3gap2_Bprime.py`, `c3gap2_verify_tip.py`, `c3gap2_heavy.py`, `c3gap2_heavy2.py`, `c3gap2_E.py`,
`c3gap2_FULLCHAIN.py` (end-to-end), `c3gap2_ties.py`, `c3gap2_ties2.py`, `c3gap2_tieproof.py`. Shared
`mn.py`. Builds on `2026-06-14-threerow-c3-Jstar-even.md`, `2026-06-14-numberlemma-general-c.md`.
