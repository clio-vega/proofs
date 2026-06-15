# The Boundary Lemma (Gap 3) for the three-row family: `c = 1, 2` closed; `c = 3` top closed

**Date:** 2026-06-15 (prove session)
**Status:** The boundary residual flagged identically in the `c=1,2,3` interior theorems is
**closed completely for `c = 1` and `c = 2`** — turning those two into *fully complete* theorems —
and **partially for `c = 3`** (the top boundary index `j=b+c` is closed for `a` even; the remaining
`c=3` boundary indices are reduced to an explicit valuation bound and verified `m ≤ 45`). A single
uniform engine — the **factor-in-product principle** — drives every closed case, so the route to
"all `c`" is now concrete.

Companion to `2026-06-13-threerow-c1-Jstar-even.md`, `2026-06-13-threerow-c2-Jstar-even.md`,
`2026-06-14-threerow-c3-Jstar-even.md`, `2026-06-15-compensation-lemma-B-proved.md`. Those proved
the **interiors**; this note closes the **boundary** (§7 / §6 residuals).

---

## 0. Setup and the object

Fix `m ≥ 1`, `λ = (a,b,c) ⊢ 2m` with `a ≥ b ≥ c ≥ 1`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`,

> `G_j := ⟨s_λ , e₂^j h₁^{2(m−j)}⟩ = C(m,j)\,M_j ∈ ℤ_{≥0}`  (`0 ≤ j ≤ m`),
> `val(j) = j + 2\,v₂(G_j)` (`= ∞` when `G_j = 0`),  `V = min_j val(j)`,  `J* = {j : val(j)=V}`.

The leading-`π` dichotomy (`π=1+i`, Lean-checked): `|J*|` odd `⟹ G_λ(i) ≠ 0`; the residual ties are
the shapes with `|J*|` even. Throughout `s₂` = binary digit sum, `v₂` = 2-adic valuation,
`v₂(k!) = k − s₂(k)` (Legendre), and **carries**`(x,y) := s₂(x)+s₂(y)−s₂(x+y) ≥ 0` is the number of
carries in the base-2 addition `x+y` (Kummer: `carries(x,y) = v₂C(x+y,x)`). Note `val(j) ≡ j
(mod 2)`.

`M_j` is supported on `0 ≤ j ≤ b+c`. The interior theorems prove `j₀ := min J*` lies in `0 ≤ j₀ ≤ b`
(in fact `j₀∈{0,3}`), and pin `J*` to an interior 2-adic box `j₀ + {0,2,…,2c}`. They each leave
**one** residual: that no **boundary index** `j ∈ {b+1,…,b+c}` is an extra minimizer.

> **Boundary Lemma (Gap 3).** Suppose `b` is large enough that the interior box lies in `[0,b]`
> (i.e. `j₀ + 2c ≤ b`; the finitely many smaller `b` are the per-family boundary-tie cases handled
> in each interior note). Then for every boundary index `j ∈ {b+1, …, b+c}`,
> `  val(j) > V.`

### 0.1 Why the naïve reduction fails (and the right one)

Since `v₂(G_j) ≥ 0`, `val(b+i) ≥ b+i ≥ b+1`. One is tempted to prove `V ≤ b`. **This is false:** e.g.
`λ=(23,8,1)`, `m=16` has `V = 21 > 8 = b`. Both `val(b+i)` *and* `V` grow with the valuations of
binomials, and they co-vary; there is no shortcut around comparing the two valuations directly.

The correct reduction is exact. By parity, write the threshold relative to `val(0)`: the interior
theorems give `V = val(0) − θ` with `θ = 0` when `j₀=0` and `θ = 3` when `j₀ = 3` (the "forced
descent" `val(3)=val(0)−3` of the `a`-odd families). So

> **`val(b+i) > V`  ⟺  `Δ(b+i) := val(b+i) − val(0) > −θ`,**

and since `val(b+i)−val(0) ≡ b+i (mod 2)`, the strict inequality sharpens to `Δ(b+i) ≥ 2−θ` when
`b+i` has the parity of `j₀`, and is automatic on the opposite parity *unless* `Δ` dips below `−θ`.
Concretely we must prove:

| family | `θ` | required bound |
|---|---|---|
| `c=1`, `a` even (`j₀=0`) | 0 | `Δ(b+1) ≥ 2` |
| `c=1`, `a` odd (`j₀=3`)  | 3 | `Δ(b+1) ≥ −2` |
| `c=2` (`j₀=0` always)    | 0 | `Δ(b+1),Δ(b+2) > 0` |
| `c=3`, `a` even (`j₀=0`) | 0 | `Δ(b+1),Δ(b+2),Δ(b+3) > 0` |
| `c=3`, `a` odd (`j₀=3`)  | 3 | `Δ(b+1),Δ(b+2),Δ(b+3) > −3` |

---

## 1. Three engines

### 1.1 The 3-row hook formula (`val(0)`)

> **Lemma H.** For `λ=(a,b,c)⊢2m`,
> `  f^λ = M_0 = \dfrac{(2m)!\,(a−c+2)(a−b+1)(b−c+1)}{(a+2)!\,(b+1)!\,c!}.`

*Proof.* Direct from the hook-length formula; the column lengths are `λ'_j = 3` (`j≤c`), `2`
(`c<j≤b`), `1` (`b<j≤a`), and the three row products telescope to `(a+2)!/[(a−c+2)(a−b+1)]`,
`(b+1)!/(b−c+1)`, `c!`. ∎ *(Verified all `λ`, `c≤5`, `m≤30`, 1662 cases.)*

### 1.2 The top boundary value (`a`-free)

> **Lemma T.** `M_{b+c} = \dfrac{(b−c+1)\,(b+c)!}{(b+1)!\,c!} = (b−c+1)\dprod_{k=2}^{c}(b+k)\big/c!.`
> Consequently, writing `P := m−b−c ≥ 0` (so `a = 2P+b+c`),
> `  Δ(b+c) = (b+c)+4 + 2 s₂(P) − 2 s₂(a+2) − 2 v₂(a−c+2) − 2 v₂(a−b+1).`

*Proof of the `Δ` formula.* `G_{b+c}/G_0 = C(m,b+c)M_{b+c}/f`. Using Lemma H and Lemma T, the
binomial `(b+c)!` and the factors `(b+1)!c!(b−c+1)` cancel, leaving
`G_{b+c}/G_0 = m!\,(a+2)! \big/ [\,P!\,(2m)!\,(a−c+2)(a−b+1)\,]`. Apply Legendre to the four
factorials; `m + (a+2) − P − 2m = 2` after `a+b+c=2m`. ∎ *(Both formulas verified `c≤5`, `m≤40`:
1520 + 2995 cases.)*

### 1.3 The factor-in-product principle (the workhorse)

This is the one elementary fact that closes every case below.

> **Lemma F.** Let `R ≥ 0`, `Q ≥ 4`. Among the `Q` consecutive integers `R+1, …, R+Q` any single
> one, say `R+t₀` (`1 ≤ t₀ ≤ Q`), satisfies
> `  v₂\Big(\dprod_{t=1}^{Q}(R+t)\Big) − v₂(R+t₀) \ \ge\ 1.`
> More generally the same holds for `Q ≥ 3` provided some index `t₁ ≠ t₀` in `[1,Q]` has `R+t₁`
> even.

*Proof.* `R+t₀` is one factor of the product, so `v₂(∏) ≥ v₂(R+t₀) + v₂(∏_{t≠t₀})`. Among `Q ≥ 4`
consecutive integers there are at least two even ones; at least one even integer `R+t₁` has
`t₁ ≠ t₀`, contributing `v₂(R+t₁) ≥ 1` to `∏_{t≠t₀}`. ∎ *(Verified `Q∈[4,80]`, `R<2000`; the bound
fails for some `R` at `Q=3`, which is exactly why each family needs its smallest `b` excluded.)*

The mechanism: a carry-sum `carries(R,Q)` plus a factorial valuation `v₂(Q!)` *is* `v₂(∏_{t=1}^Q
(R+t))` (Kummer), and the dangerous subtracted term `−2v₂(\text{linear in } a)` always equals
`−2v₂(R+t₀)` for a factor `R+t₀` of that product. Lemma F then yields a clean surplus.

---

## 2. `c = 1` — COMPLETE

`λ=(a,b,1)`, `a+b = 2m−1` (so `a,b` opposite parity). Residual range: `b ∉ {1,2,4}` (i.e. `b≥3`
for `a` even with `b` odd, `b ≥ 6` for `a` odd with `b` even); the excluded `b` are the
boundary-tie families of the `c=1` note §6.

Only `j = b+1` is a boundary index (`M_{b+1} = b`, `M_j = 0` for `j ≥ b+2`).

> **Proposition 1 (clean `Δ`).** `Δ(b+1) = (b−3) + 2 v₂(a+2) + 2\big(s₂(m−b) − s₂(a)\big).`

*Proof.* `Lemma H` for `c=1` is `f = (2m)(2m−1)\,C(2m−2,a)(a−b+1)/[(a+2)(b+1)]`, giving
`val(0) = 2 + 2v₂(m) + 2v₂C(2m−2,a) + 2v₂(a−b+1) − 2v₂(a+2) − 2v₂(b+1)`. With `M_{b+1}=b`,
`val(b+1) = (b+1) + 2v₂(b) + 2v₂C(m,b+1)`. Subtracting and using (i) `b·C(m,b)=m·C(m−1,b−1)`,
(ii) the symmetry `C(2(m−1),a) = C(2(m−1),b−1)` (as `a = 2(m−1)−(b−1)`), (iii) `a−b+1 = 2(m−b)` so
`v₂(a−b+1) = 1+v₂(m−b)`, collapses everything to two binomials of equal lower index:
`Δ(b+1) = (b−3) + 2v₂(a+2) + 2\big[v₂C(m−1,b−1) − v₂C(2(m−1),b−1)\big]`. Kummer
(`v₂C(M,k)−v₂C(2M,k) = s₂(M−k) − s₂(2M−k)`) with `M=m−1`, `k=b−1`, `M−k=m−b`, `2M−k=a`, finishes. ∎
*(Verified all `λ=(a,b,1)`, `m ≤ 60`, 1711 cases, 0 mismatch.)*

Now `a = 2(m−b) + (b−1)` gives `s₂(m−b) − s₂(a) = γ − s₂(b−1)` with `γ := carries(2(m−b),\,b−1) ≥ 0`,
so

> `Δ(b+1) = (b−3) + 2 v₂(a+2) + 2γ − 2 s₂(b−1)`.

Recall **Lemma A** (digit-sum envelope): `2s₂(k) ≤ k−1` for `k ≥ 4`; `=k+1` only at `k∈{1,3}`.

- **`a` odd** (`b` even `≥ 6`, need `Δ ≥ −2`; `Δ` is odd `⟹` need `≥ −1`): `v₂(a+2)=0`,
  `b−1` odd `≥ 5`, so `2s₂(b−1) ≤ b−2`. Then `Δ ≥ (b−3) + 2γ − (b−2) = 2γ − 1 ≥ −1`. ✓
- **`a` even** (`b` odd `≥ 3`, need `Δ ≥ 2`; `Δ` is even): `v₂(a+2) ≥ 1`.
  - `b ≥ 5`: `b−1` even `≥ 4`, `2s₂(b−1) ≤ b−2`, so `Δ ≥ (b−3)+2+2γ−(b−2) = 1+2γ ≥ 1`; even `⟹ ≥ 2`. ✓
  - `b = 3`: `Δ = 2v₂(a+2) + 2γ − 2`. Here `a = 2m−4`, `v₂(a+2) = 1 + v₂(m−1)`. If `m` is odd,
    `v₂(a+2) ≥ 2 ⟹ Δ ≥ 2`. If `m` is even, `v₂(a+2)=1` and `γ = carries(2(m−3),2) ≥ 1` (the bit-1
    of `2(m−3)` is set because `m−3` is odd), so `Δ = 2γ ≥ 2`. ✓

Hence `val(b+1) > V` for all `b ∉ {1,2,4}`: **the `c=1` boundary lemma is proven, and the three-row
`c=1` theorem is now complete.** *(Conclusion re-verified directly, `m ≤ 50`, 1036 cases, margin ≥ 2.)*

---

## 3. `c = 2` — COMPLETE

`λ=(a,b,2)`, `a+b = 2m−2` (so `a,b` same parity), `j₀ = 0`, `θ=0`. Residual `b ≥ 3` (`b=2` is the
note's §6 family). Set `P := m−b−2`. Boundary: `M_{b+1} = (b−1)W₁/2` with
`W₁ = a(b+2)−b(b+1) = 2(b+2)(P+1) + b`, and `M_{b+2} = (b−1)(b+2)/2`. Note `a−b+1` is **odd** here,
so `v₂(a−b+1) = 0` throughout.

### 3.1 Top index `j = b+2`

By Lemma T, `Δ(b+2) = (b+6) + 2 s₂(P) − 2 s₂(a+2) − 2 v₂(a)` (`a−c+2=a`). With
`a+2 = 2P + (b+4)`, Kummer gives `s₂(a+2) = s₂(P)+s₂(b+4)−δ`, `δ := carries(2P,\,b+4)`, hence

> `Δ(b+2) = (b+4) − 2 s₂(b+4) + 2δ − 2 v₂(a) \;+\; 2`.

- **`a,b` odd:** `v₂(a)=0`; `b+4` odd `≥ 7`, so `2s₂(b+4) ≤ b+3 < b+4`, giving `Δ(b+2) ≥ 2`. ✓
- **`a,b` even** (`b ≥ 4`): put `Q := b/2 + 2 ≥ 4`. Then `b+4 = 2Q`, `δ = carries(2P,2Q) =
  carries(P,Q)`, `v₂(a) = 1 + v₂(P+Q−1)`, and `(b+4)−2s₂(b+4) = 2(Q−s₂(Q)) = 2 v₂(Q!)`. So
  `Δ(b+2) = 2\big[v₂(Q!) + carries(P,Q) − v₂(P+Q−1)\big] = 2\big[v₂(∏_{t=1}^{Q}(P+t)) −
  v₂(P+Q−1)\big]`. Since `P+Q−1` is the `t=Q−1` factor, **Lemma F** (`Q≥4`) gives `Δ(b+2) ≥ 2`. ✓
  *(Identity verified `m ≤ 45`, 400 cases.)*

### 3.2 Subtop index `j = b+1`

By the same hook+Kummer reduction (`a−b+1` odd),
`Δ(b+1) = (b+3) + 2 s₂(P+1) − 2 s₂(a+2) + 2 v₂(W₁) − 2 v₂(a)`, and with `a+2 = 2(P+1)+(b+2)`,
`s₂(a+2) = s₂(P+1)+s₂(b+2)−δ₁`, `δ₁ := carries(2(P+1),\,b+2)`:

> `Δ(b+1) = (b+3) − 2 s₂(b+2) + 2δ₁ + 2 v₂(W₁) − 2 v₂(a)`.

- **`a,b` odd:** `v₂(W₁)=v₂(a)=0`; `b+2` odd `≥ 5`, `2s₂(b+2) ≤ b+1`, so `Δ(b+1) ≥ 2+2δ₁ ≥ 2`. ✓
- **`a,b` even** (`b = 2β`, `β ≥ 2`, need `Δ(b+1) ≥ 1`, parity odd): write `X := P+1`. Then
  `a = 2(X+β)`, `W₁ = 2[(b+2)X + β]`, `δ₁ = carries(X, β+1)`, `b+3−2s₂(b+2) = 2β+3 − 2s₂(β+1)`.
  Collecting and using `2β+3−2s₂(β+1) + 2 carries(X,β+1) = 1 + 2 v₂(∏_{t=1}^{β+1}(X+t))` (Kummer,
  `(β+1)−s₂(β+1)=v₂((β+1)!)`):

  > `Δ(b+1) = 1 + 2\Big[\,v₂\big(∏_{t=1}^{β+1}(X+t)\big) + v₂(2(β+1)X+β) − v₂(X+β)\,\Big].`

  The bracket is `≥ 0` because **`X+β` is the `t=β` factor** of `∏_{t=1}^{β+1}(X+t)` (so the first
  minus the third term is `≥0`) and `v₂(2(β+1)X+β) ≥ 0`. Hence `Δ(b+1) ≥ 1`. ✓ *(Identity +
  certificate verified `m ≤ 45`, 441 cases.)*

So `val(b+1),val(b+2) > val(0) = V` for all `b ≥ 3`: **the `c=2` boundary lemma is proven, and the
three-row `c=2` theorem is now complete.** *(All boundary indices, all parities, `m ≤ 45`, 1764
cases, 0 violation.)*

---

## 4. `c = 3` — top closed for `a` even; remainder reduced & verified

`λ=(a,b,3)`, `a+b = 2m−3` (`a,b` opposite parity), `j₀ ∈ {0,3}` (`a` even `→ 0`, `a` odd `→ 3`),
`P := m−b−3`. Box interior requires `b ≥ 6` (`a` even) / `b ≥ 9` (`a` odd); smaller `b` are the
`c=3` note's boundary-tie families.

### 4.1 Top index `j = b+3`, `a` even — COMPLETE

`a−c+2 = a−1` is **odd** (`a` even), so `v₂(a−1)=0`; and `a−b+1 = 2(P+2)`. Lemma T plus the carry
expansion `a+2 = 2P+(b+5)` give

> `Δ(b+3) = (b+5) − 2 s₂(b+5) + 2Γ − 2 v₂(P+2)`,  `Γ := carries(2P,\,b+5)`.

`b` is odd (`b ≥ 7`), so `b+5` is even; put `Q' := (b+5)/2 ≥ 6`. Then `b+5 = 2Q'`, `Γ =
carries(P,Q')`, `(b+5)−2s₂(b+5) = 2 v₂(Q'!)`, hence

> `Δ(b+3) = 2\Big[\,v₂\big(∏_{t=1}^{Q'}(P+t)\big) − v₂(P+2)\,\Big].`

`P+2` is the `t=2` factor of the product, and `Q' ≥ 6`, so **Lemma F** gives `Δ(b+3) ≥ 2 > 0`. ✓
*(Identity verified `m ≤ 45`, 342 cases; margin `≥ 2`.)*

This is the first confirmation that the factor-in-product engine handles the genuinely harder `c=3`
boundary — the route to all `c` for the top index is exactly this (Lemma T is uniform in `c`; only
the parity of `a−c+2`, `a−b+1` shifts which factor `P+t₀` absorbs the deficit).

### 4.2 The remaining `c=3` indices — precise residual

For the subtops `j=b+1,b+2` the closed forms carry an extra polynomial factor in `a`:
`M_{b+2} = (b−2)(b+2)\,N₂/6` with `N₂ = a(b+3) − (b²+2b+3)`, and `M_{b+1} = (b−2)(a−b+1)\,N₁/12`
with `N₁ = ab²+5ab+6a − b³−b²−4b−6`. Telescoping from the top (`G_{j}/G_{j+1} =
\tfrac{j+1}{m−j}\tfrac{M_{j}}{M_{j+1}}`) gives the exact

> `Δ(b+2) = Δ(b+3) − 1 + 2 v₂(N₂) − 2 v₂(P+1)`,
> `Δ(b+1) = Δ(b+2) − 3 + 2 v₂(a−b+1) + 2 v₂(N₁) − 2 v₂(P+2) − 2 v₂(N₂)`.

These reduce the lemma to bounding `v₂(N₁), v₂(N₂)` (valuations of explicit quadratic/cubic forms in
`a`) against `v₂(P+1),v₂(P+2)`. The factor-in-product engine does not yet apply verbatim because
`N₁,N₂` are not single factors of a consecutive-integer product. **This is the precise residual.**

**Status (honest).** `Δ(b+i) > −θ` for all `c=3` boundary indices, both parities, holds with margin
`≥ 2` over the box-interior range: verified `m ≤ 39` (1368 cases, 0 violation) here, and `m ≤ 79`
in the `c=3` interior note. No boundary tie occurs (strict), confirming the interior theorem needs
no correction. What remains for a *hand* proof of `c=3`: a bound `v₂(N₂) ≥ v₂(P+1) − 1` (would close
`j=b+2`, both parities) and the analogous bound for `N₁`; and the `a`-odd top, where **both**
`v₂(a−1)` and `v₂(a−b+1)` are positive and two factors of the product must absorb the two deficits
(a two-factor Lemma F).

---

## 5. Verification summary

All scripts under `projects/code/threerow-boundary/` (reuse `threerow-c1/mn.py`).

| claim | range | cases | result |
|---|---|---|---|
| Lemma H (hook formula) | `c≤5`, `m≤30` | 1662 | 0 mismatch |
| Lemma T (`M_{b+c}` + `Δ(b+c)`) | `c≤5`, `m≤40` | 1520+2995 | 0 mismatch |
| Lemma F (factor-in-product) | `Q≤80`, `R<2000` | — | 0 fail (`Q≥4`) |
| `c=1` Prop 1 + boundary lemma | `m≤60` / `m≤50` | 1711 / 1036 | 0 / margin ≥2 |
| `c=2` top + subtop identities | `m≤45` | 400+441 | 0 mismatch |
| `c=2` full boundary `Δ>0` | `m≤45`, all parities | 1764 | 0 violation |
| `c=3` `a`-even top identity | `m≤45` | 342 | 0, margin ≥2 |
| `c=3` full boundary (box interior) | `m≤39` | 1368 | 0 violation |

---

## 6. What this closes, and the path forward

- **`c = 1` and `c = 2` three-row `d=4` families are now COMPLETE theorems** — interior + boundary,
  no machine-checked residual remains. `J*` is exactly the interior 2-adic box; `|J*|` is even on
  every tie.
- The **Boundary Lemma reduces** (exactly, not heuristically) to explicit elementary inequalities
  `Δ(b+i) ≥ 2−θ`, all of the shape *"positive linear in `b` + carries − (deficit) ≥ 0"*, where the
  deficit is `v₂` of one or two factors of a consecutive-integer product.
- The **uniform engine is Lemma F** (factor-in-product): the subtracted `v₂(\text{linear in }a)`
  is always a factor `P+t₀` of the consecutive product `∏(P+t)` that the carry-sum builds, and
  `≥4` consecutive integers always carry a spare factor of 2. This closes `c=1` (both parities),
  `c=2` (both indices, both parities), and the `c=3` top (`a` even) **with the same one-line bound**.
- **Remaining for the full win (`c=3`):** bound `v₂(N₁),v₂(N₂)` for the two subtop indices and run a
  two-factor Lemma F for the `a`-odd top. Both are now isolated, finite-shape, and verified — a
  bounded follow-on (and the natural LEAN target alongside the now-complete `c=1,2` families).

### Files
`projects/code/threerow-boundary/`: `boundary_forms.py` (Lemmas T, H + boundary `M` fits),
`c1_clean.py` (`c=1` Prop 1), `c2_final.py` (`c=2` identities + full check), `c3_aeventop.py`
(`c=3` top), `lemma_only.py` (Lemma F). Shared `mn.py` (Murnaghan–Nakayama).
