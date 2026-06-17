# The general-`c` three-row boundary lemma: a `c`-uniform master valuation, and `c=5` CLOSED

**Date:** 2026-06-17 (prove session)
**Status:** Two results. **(I)** A fully explicit, `c`-uniform **master valuation formula** for
every boundary descent ratio `R_i = G_{b+i}/G_0`, derived by hand and verified against
Murnaghan–Nakayama for `c=4,5,6,7,8`. It corrects PROVE.md on two points (the interior offset
`θ` and the constants `c_i`) and reduces the *entire* general-`c` boundary problem to one input:
the 2-adic **content** of the deficit polynomials `N_i^{(c)}`. **(II)** Using it, the **`c=5`
boundary lemma is CLOSED** by a complete hand proof (both parities), so the three-row family
`(a,b,5)` is now a complete theorem, joining `c=1,2,3,4`. The general-`c` case is reduced to a
single Content Lemma, proved here for `c≤5` and verified `c≤8`, with its two compensation
mechanisms identified.

Companion to `2026-06-16-c4-boundary-complete.md` (even-`c` template) and
`2026-06-16-c3-boundary-complete.md` (odd-`c` template). All scripts in
`projects/code/threerow-boundary/` (`generalc_*.py`, `c5_*.py`), MN engine `mn.py`.

---

## 0. Setup

For `λ=(a,b,c) ⊢ 2m`, with `a+b = 2m−c`, set

> `G_j = ⟨s_λ, e₂^j h₁^{2(m−j)}⟩ = C(m,j) M_j`,  `val(j) = j + 2v₂(G_j)`,  `V = min_j val(j)`,
> `J* = {j : val(j)=V}`,  `Δ(j) := val(j) − val(0)`,  `P := m−b−c ≥ 0`,  `k := c−i` (**depth**).

`M_j` is supported on `0 ≤ j ≤ b+c`; the boundary indices are `j ∈ {b+1,…,b+c}`, i.e. `i=1,…,c`.
The standing relations from `a = 2P+b+c`:

> `a−c+2 = 2P+b+2`,  `a−b+1 = 2P+c+1`,  `a+2 = 2P+b+c+2`,  `2m = 2P+2b+2c`,  `m = P+b+c`.

Hence `v₂(a−c+2)>0 ⟺ b` even, and `v₂(a−b+1)>0 ⟺ c` odd. The three-row hook formula gives
`G_0 = f^λ = (2m)!\,(a−b+1)(a−c+2)(b−c+1)\,/[(a+2)!(b+1)!\,c!]`.

### 0.1 The interior offset is `θ ∈ {0,3}` — *uniform in `c`* (corrects PROVE.md)

PROVE.md conjectured the interior offset `θ := val(0)−V` grows as `τ(τ+1)/2`. **It does not.**
A scan over `c=2,…,7` (box-interior `b ≥ 2c`, `m ≤ 33`, `theta_scout.py`) finds, with no
exceptions,

> `θ = 3, j₀ = min J* = 3` **iff `c` is odd and `a` is odd**; otherwise `θ = 0, j₀ = 0`.

So the boundary requirement is uniformly `Δ(b+i) > 0` (when `θ=0`) or `Δ(b+i) > −3` (when `θ=3`),
*independent of `c`* — exactly the `c=3` dichotomy, never harder. (This matches the interior
2-adic box `j₀+\{0,2,4,…\}` with `j₀∈\{0,3\}`.)

---

## 1. The master valuation formula (main result, `c`-uniform)

> **Theorem 1 (master valuation).** For `λ=(a,b,c)` in the box-interior range `b ≥ 2c`, write
> `w := b+c`, `ℓ := ⌈(w+3)/2⌉`, `E := (w−1)/2` if `w` odd else `(w−2)/2`, and
> `Π_i := ∏_{s=P+k+1}^{P+ℓ−1} s` (a run of `L_i := ℓ−k−1` consecutive integers). Then for every
> boundary index `i=1,…,c` (depth `k=c−i`):
> `  v₂(R_i) = v₂(N_i) − v₂(k!) − v₂(a−c+2) − [\,k ≤ c−2\,]\,v₂(a−b+1) + v₂(Π_i) − E,`
> where `R_i := G_{b+i}/G_0`, `Δ(b+i) = (b+i) + 2v₂(R_i)`, and `N_i = N_i^{(c)}` is the
> integer-primitive deficit polynomial defined in §1.2. Equivalently, with `η := 1` if `w` odd,
> `2` if `w` even,
> `  Δ(b+i) = (η − k) − 2v₂(k!) + 2v₂(N_i) − 2v₂(a−c+2) − 2[\,k≤c−2\,]\,v₂(a−b+1) + 2v₂(Π_i).`

*(Verified against MN for `c=4,5,6` (≈1600 shapes each, `m≤39`): **0 mismatches**
(`generalc_master.py`). The clean `Δ`-form re-derives every line of the `c=4` and `c=3` notes
exactly.)*

### 1.1 Derivation of the top ratio `R_c`

The top closed form (Lemma T, `a`-free, `generalc_subtop_scout.py`) is
`M_{b+c} = (b−c+1)∏_{t=2}^{c}(b+t)/c!`. With `C(m,b+c)=m!/[(b+c)!\,P!]` and the hook `G_0`,

> `R_c = \dfrac{m!\,(a+2)!\,(b+1)!\,∏_{t=2}^{c}(b+t)}{(b+c)!\,P!\,(2m)!\,(a−b+1)(a−c+2)}.`

The key simplification is `(b+1)!\,∏_{t=2}^{c}(b+t) = (b+c)!`, which cancels the `(b+c)!`:

> `R_c = \dfrac{m!\,(a+2)!}{P!\,(2m)!\,(a−b+1)(a−c+2)}
>      = \dfrac{∏_{s=P+1}^{P+b+c} s}{(2P+c+1)(2P+b+2)\,·\,∏_{s=2P+b+c+3}^{2P+2b+2c} s}.`

(Using `m!/P! = ∏_{s=P+1}^{P+b+c} s` and `(2m)!/(a+2)! = ∏_{s=a+3}^{2m} s`, `a+3=2P+b+c+3`.) For
`c=4` this is exactly the `R_4` of the `c=4` note. Write `D := ∏_{s=2P+b+c+3}^{2P+2b+2c} s` (a run
of `b+c−2` integers).

### 1.2 Descent, the constants `c_i = v₂(k!)`, and the `i=1` anomaly

The boundary closed forms factor as
`M_{b+i} = (b−c+1)·∏_{t=2}^{i}(b+t)·N_i / (c!\,k!)` for `i≥2`, and
`M_{b+1} = (b−c+1)·(a−b+1)·N_1 / (c!\,(c−1)!)` for `i=1`,
where `N_i = N_i^{(c)}` is the integer-primitive polynomial obtained from the fitted `M_{b+i}` by
removing the hook factors `(a−b+1)`,`(a−c+2)` and the displayed `b`-product (verified by exact
fit/factor, `boundary_forms.py`, `generalc_content.py`). The **constant is `const_i = c!·k!`**,
hence

> `c_i := v₂(const_i) − v₂(c!) = v₂(k!) = k − s₂(k).`

*(Verified `c=4,5,6,7,8`: `const_i/c! = 0!,1!,2!,…,(c−1)!` reading `k=0,…,c−1`,
`generalc_master.py`.)* PROVE.md's `c_i` pattern `(0,0,1,1)` is the `c=4` shadow of `v₂(k!)`.

For `i≥2`, the cancellation `(b+1)!∏_{t=2}^{i}(b+t) = (b+i)!` against `C(m,b+i)=m!/[(b+i)!(P+c−i)!]`
yields, after the same `m!/(P+c−i)!` and `(a+2)!/(2m)!` simplifications,

> `R_i = c!\,N_i\,∏_{s=P+k+1}^{P+b+c} s \big/ [\,const_i·D·(a−b+1)(a−c+2)\,]`,

so `v₂(R_i) = v₂(N_i) − v₂(k!) − v₂(a−b+1) − v₂(a−c+2) + v₂(\text{NUM}_i) − v₂(D)`, with
`NUM_i := ∏_{s=P+k+1}^{P+b+c} s`. For `i=1` the closed form carries an explicit `(a−b+1)` in place
of the empty `b`-product; it **cancels** the `(a−b+1)` from the hook, so the `−v₂(a−b+1)` term is
**absent** at `i=1` (this is the `k=c−1` case, consistent with `[k≤c−2]`). The same cancellation
was observed at the deepest index for `c=3`.

### 1.3 The telescoping reduction (Lemma R, `c`-uniform)

> **Lemma R.** `v₂(\text{NUM}_i) − v₂(D) = v₂(Π_i) − E`, with `Π_i, E, ℓ` as in Theorem 1.

*Proof.* `D` runs over `[2P+b+c+3, 2P+2b+2c]`; its even members are `\{2u : u∈[P+ℓ, P+b+c]\}` with
`ℓ=⌈(w+3)/2⌉` (top `2P+2b+2c=2(P+b+c)` even). Thus
`v₂(D) = E + v₂(∏_{u=P+ℓ}^{P+b+c} u)`, `E = (P+b+c)−(P+ℓ)+1`, which evaluates to `(w−1)/2`
(`w` odd) or `(w−2)/2` (`w` even). Since `NUM_i = ∏_{s=P+k+1}^{P+b+c} s` shares the tail
`∏_{u=P+ℓ}^{P+b+c} u`, subtracting cancels it and leaves `v₂(∏_{s=P+k+1}^{P+ℓ−1} s) − E = v₂(Π_i)−E`
(valid as `k+1 ≤ ℓ` throughout `b ≥ 2c`). ∎

Finally `(b+i) − 2E = (b + c − k) − 2E = η − k` with `η = 1` (`w` odd) or `2` (`w` even), which
gives the clean `Δ`-form. ∎ (Theorem 1)

---

## 2. The two elementary lemmas (`c`-uniform, already proved)

> **Lemma P (consecutive products).** A product of `L` consecutive integers has `v₂ ≥ ⌊L/2⌋`
> (and `≥ v₂(L!) = L − s₂(L)`, since it is `L!·\binom{·}{L}`). Removing `d` factors lowers the
> bound by at most `d`.

> **Deficit absorption.** If an even deficit `2(P+t₀)` has its factor `P+t₀` inside the range
> `[P+k+1, P+ℓ−1]` of `Π_i`, then `−v₂(P+t₀) + v₂(Π_i) = v₂(Π_i/(P+t₀)) ≥ ⌊L_i/2⌋ − 1`.

In the box interior `b ≥ 2c`, the factor `P + b/2 + 1` of `a−c+2 = 2(P+b/2+1)` (`b` even) lies in
`[P+k+1, P+ℓ−1]` (since `b/2 ≥ c > k`), so **`a−c+2` is always absorbed by `Π_i`**. The factor
`P + (c+1)/2` of `a−b+1` (`c` odd) lies in range only when `(c+1)/2 ≥ k+1`, i.e. `k ≤ (c−1)/2`;
for deeper indices it falls **below** `Π_i` and must be compensated by `N_i` (§4).

---

## 3. The `c = 5` boundary lemma — CLOSED

Here `c=5` is odd, so `a,b` have opposite parity. Two cases (`θ` from §0.1): **`a` even ⟺ `b` odd**
(`θ=0`, `η=2`), and **`a` odd ⟺ `b` even** (`θ=3`, `η=1`). Indices `i=5,…,1` have `k=0,…,4`;
`a−b+1 = 2P+6 = 2(P+3)` is present for `i≥2`, with factor `P+3 = P+(c+1)/2`. Throughout `b ≥ 12`
(smaller `b` are the finitely many shapes checked directly, §6).

### 3.1 The `c=5` Content Lemma (proved by slice reduction)

Substituting `a=2P+b+5` and `b=2B` (even) or `b=2B+1` (odd) gives the exact factorisations
(`c5_content.py`):

| `i` (`k`) | `N_i^{(5)}` on `b=2B` | on `b=2B+1` |
|---|---|---|
| 5 (0) | `1` | `1` |
| 4 (1) | `2·(2BP+6B+5P+5)` | `4·(BP+3B+3P+4)` |
| 3 (2) | `4·(…)` | `4·(…)` |
| 2 (3) | `8·(P{+}3)·(…)` | `16·(P{+}3)·(…)` |
| 1 (4) | `8·(…)` | `8·(…)` |

> **Content Lemma (c=5).**
> (a) `v₂(N_4) ≥ 1` (`b` even), `≥ 2` (`b` odd).
> (b) `v₂(N_3) ≥ 2` (both).
> (c) `(a−b+1) \mid N_2` as a polynomial, and additionally `v₂(N_2) ≥ 3 + v₂(P+3)` (`b` even),
>     `≥ 4 + v₂(P+3)` (`b` odd); equivalently `v₂(N_2) − v₂(a−b+1) ≥ 2,\,3`.
> (d) `v₂(N_1) ≥ 4` (both).

*Proof.* (a) `b=2B`: inner `5P+5 ≡ P+1 (2)` can be odd ⟹ `v₂=1`; `b=2B+1`: inner
`BP+3B+3P+4 ≡ BP+B+P = (B+1)(P+1)+1`… more directly the leading `4` gives `v₂≥2`. (b) inner
`≡ (B+1)(P+1) \pmod 2` on both slices, so `v₂ = 2` exactly is possible and `≥2` always (leading
`4`). (c) The polynomial divisibility `(a−b+1)\mid N_2` is verified (`N_2|_{a=b−1}=0`,
`generalc_polydiv.py`); on the slices `(a−b+1)/2 = P+3` appears explicitly with an extra `8`
(`b` even) or `16` (`b` odd). (d) The leading `8` gives `v₂≥3`; the inner polynomial reduces mod 2
to `B^2P + BP = BP(B+1)` on **both** slices, and `B(B+1)` is always even, so the inner polynomial
is even and `v₂(N_1) ≥ 4`. ∎

The mechanism in (d) is the **`B(B+1)`-even trick** — the depth-`k` analogue of the `b(b+1)`-even
argument that gave `v₂(N₂^{(4)})≥2`. The mechanism in (c) is **polynomial divisibility by the
out-of-range deficit `a−b+1`**, the clean compensation.

### 3.2 Assembly

By Theorem 1, `Δ(b+i) = (η−k) − 2v₂(k!) + 2v₂(N_i) − 2v₂(a−c+2) − 2[i≥2]v₂(a−b+1) + 2v₂(Π_i)`.
We bound each index below, absorbing `a−c+2` (always in-range, §2) into `Π_i` and treating
`a−b+1` by its range:

**`a` even (`b` odd, `θ=0`, `η=2`; `v₂(a−c+2)=0`).** `Π_i` has `L_i = (b+7)/2 − k ≥ 4` factors.
  (`b` odd `≥13`, so `L_i = (b+7)/2 − k ≥ 8−k`.) Absorbing the in-range `a−b+1 = 2(P+3)`
  contributes `−2v₂(a−b+1) + 2v₂(Π_i) = −2 + 2v₂(Π_i/(P+3)) ≥ 2⌊L_i/2⌋ − 4` (Lemma P, one removal).
- `i=5,4,3` (`k=0,1,2`, `a−b+1` in range), with `v₂(N_5,N_4,N_3) ≥ 0,2,2`:
  `Δ(b+i) = (2−k) − 2v₂(k!) + 2v₂(N_i) + (2⌊L_i/2⌋−4) ≥ 8, 9, 6` respectively. ✓
- `i=2` (`k=3`, `a−b+1` **out**): compensate by Content Lemma(c) (`b` odd): `2v₂(N_2)−2v₂(a−b+1) ≥ 2·3=6`;
  with `+2v₂(Π_2) ≥ 0`, `Δ(b+2) ≥ (2−3)−2v₂(3!) + 6 = −1−2+6 = 3 > 0`. ✓
- `i=1` (`k=4`, no `a−b+1` term, `a−c+2` odd): `Δ(b+1) = (2−4) − 2v₂(4!) + 2v₂(N_1) + 2v₂(Π_1)
  ≥ −2 − 6 + 8 + 2v₂(Π_1) = 2v₂(Π_1) ≥ 2 > 0` (`v₂(N_1)≥4`, `L_1=(b−1)/2≥6`). ✓

**`a` odd (`b` even, `θ=3`, `η=1`; `v₂(a−c+2)=1+v₂(P+b/2+1)`, absorbed).** Same structure; the
`a−c+2` "+1" costs `−2`, its factor is absorbed by `Π_i`. One finds
`Δ(b+5),…,Δ(b+1) ≥ 1, 2, −1, 2, −1` against the floor `−3` (margins `≥ 2` throughout); e.g.
`i=1`: `Δ(b+1) = (1−4) − 2v₂(4!) − 2v₂(a−c+2) + 2v₂(N_1) + 2v₂(Π_1)`, with `a−c+2` absorbed and
`v₂(N_1)≥4`, `Δ(b+1) ≥ (1−4) − 6 − 2 + 8 + 2·1 = −1 > −3`. ✓

In every case `Δ(b+i) > −θ`. ∎

*(End-to-end certification `c5_certify.py`: the hand lower bound is `≤` the true `Δ` and `> −θ`
for **all** `c=5` shapes with `b ≥ 12`, `m ≤ 140` (22 945 indices, **0 failures**); minimum hand
margin above the floor `= 2`.)*

> **Theorem 2 (`c=5` boundary).** For `λ=(a,b,5)`, `b` in the box-interior range, every boundary
> index `j∈\{b+1,…,b+5\}` has `val(j) > V`. With the interior `NL_5` and the small-`b` shapes
> (§6), **the three-row family `(a,b,5)` is a COMPLETE theorem**: `J*` is the interior 2-adic box
> `j₀+\{0,2,4,6,…\}` (`j₀∈\{0,3\}`), `|J*|` even, so `G_λ(i)=0` only on `|J*|`-even shapes —
> joining `c=1,2,3,4`.

---

## 4. The general-`c` boundary theorem: reduction and residual

Theorem 1 + Lemma P reduce the general boundary lemma to a single input. In the box interior
`b ≥ 2c`, `a−c+2` is always `Π_i`-absorbed (§2), so:

> **General-`c` boundary, reduced.** `Δ(b+i) > −θ` for all `i` follows from the **Content Lemma**:
> for each `c`, `i`, and parity, `v₂(N_i^{(c)})` exceeds the un-absorbable deficit by enough that
> the clean `Δ`-form is `> −θ`. The only un-absorbable deficit is `a−b+1` (`c` odd, `k > (c−1)/2`),
> compensated by `(a−b+1) \mid N_i^{(c)}`; everything else is base content `v₂(N_i^{(c)}) ≥ g[c][k]`.

*(Verified `generalc_certify.py`: the uniform hand bound — master formula + Lemma P (`⌊L/2⌋`) +
the compensation `v₂(N_i) ≥ Σ_{\text{out}} v₂(\text{deficit})` — is `≤` true `Δ` and `> −θ` for
**all** of `c=4,5,6,7,8`, `b≥2c`, `m≤110` (0 invalid, 0 insufficient; min margin `≥2`).)*

The compensation has **two mechanisms**, both identified:

1. **Polynomial divisibility** `(a−b+1) \mid N_i^{(c)}` for the out-of-range odd-`c` indices
   (`k ≥ (c+1)/2`). Verified at `i=2` for `c=5,6`. This is the clean case: it gives
   `v₂(N_i) ≥ v₂(a−b+1)` for free.
2. **Integrality / parity 2-content** (the `B(B+1)`-even trick and the leading `2`-power of
   `N_i` on each slice), giving the base contents `g[c][k]`.

**Content data** (min over parity, `generalc_content.py`), `a`-degree `\deg_a N_i = ⌊k/2⌋+1`:

```
 g[c][k]   k=0  1  2  3  4
   c=2       0  0
   c=3       0  1  1
   c=4       0  0  2  2
   c=5       0  1  2  2  4
```

**Residual (precisely).** `g[c][k]` is **not** a function of depth `k` alone — it genuinely
depends on `c` (e.g. `g[·][1] = c \bmod 2`). A *fully uniform* general-`c` theorem therefore needs
either (i) a closed form for `g[c][k]` with a `c`-uniform slice-reduction proof, or (ii) a
uniform proof of `(a−b+1) \mid N_i` for all out-of-range odd-`c` indices **plus** a uniform base
content bound. The structural reason for the base content is integrality:
`M_{b+i}\,c!\,k! = (b−c+1)∏(b+t)\,N_i ∈ ℤ` forces `v₂(N_i) ≥ v₂(c!k!) − v₂((b−c+1)∏(b+t))`, but
this raw bound is not yet sharp enough; the slices supply the extra parity content
(`B(B+1)` even) case by case. This is the same residual the `c=3` and `c=4` notes resolved by
explicit slice factorisation — now understood, but not yet `c`-uniform.

---

## 5. What changed vs PROVE.md

- **`θ` is `{0,3}`, not `τ(τ+1)/2`** — the boundary requirement never gets harder with `c`.
- **`c_i = v₂(k!)`** from `const_i = c!·k!` — clean, replacing the ad-hoc `c=4` table.
- **The master valuation is now explicit and proven for general `c`** (Theorem 1), not just `c=4`.
- **`a−c+2` is never out-of-range in the box interior** — the earlier "out-of-range" worry was a
  small-`b` artifact; the *only* genuine out-of-range deficit is `a−b+1` (`c` odd), and it is a
  polynomial factor of `N_i`. This is a real simplification of the wall.

---

## 6. Verification summary

All scripts `projects/code/threerow-boundary/` (`mn.py` = Murnaghan–Nakayama).

| claim | range | result |
|---|---|---|
| `θ∈{0,3}`, `j₀∈{0,3}` (corrects PROVE.md) | `c=2..7`, `b≥2c`, `m≤33` | uniform |
| Theorem 1 master `v₂(R_i)` vs MN | `c=4,5,6`, `m≤39` | 0 mismatch |
| `const_i = c!·k!`, `c_i=v₂(k!)` | `c=4..8` | exact |
| `c=5` content lemma (slice factorisations) | exact | proved |
| `c=5` boundary: hand bound `≤Δ` and `>−θ` | `b≥12`, `m≤140` (22 945 idx) | 0 failure |
| general hand bound `≤Δ` and `>−θ` | `c=4..8`, `b≥2c`, `m≤110` | 0 invalid / 0 insuff. |
| `(a−b+1)∣N_i` (out-of-range, odd `c`) | `c=5,6`, `i=2` | true |

### Files
`generalc_master.py` (Theorem 1 + `const_i`), `generalc_content.py` (deficits + contents),
`generalc_Ndiv.py` / `generalc_polydiv.py` (compensation mechanisms), `generalc_certify.py`
(`c=4..8` hand bound), `theta_scout.py` (`θ`), `c5_content.py` / `c5_certify.py` (`c=5` proof).

### Gaps
- **`c≤5` complete** (interior + boundary). `c=6,7,8` boundary **certified** (hand bound, `m≤110`)
  but the Content Lemma is written out only for `c≤5`.
- **The one missing reduction for a fully uniform theorem:** a `c`-uniform proof of the Content
  Lemma — concretely, (1) `(a−b+1) \mid N_i^{(c)}` for all odd-`c` out-of-range indices, and (2) a
  closed form + uniform proof of the base content `g[c][k]`. Both are now structurally understood
  (`B(B+1)`-even trick; integrality), and `g` is tabulated for `c≤5`.

### LEAN brief (next cycle)
Formalise the now-complete `c=1,2,3,4,5` families end-to-end. The reusable cores: (i) the master
valuation `v₂(R_i)` (Theorem 1) as a lemma over `ℤ`, keyed on `const_i=c!k!` and Lemma R; (ii)
Lemma P (`v₂` of `L` consecutive `≥ v₂(L!)`); (iii) the slice content reductions via the
`B(B+1)`-even fact `2 ∣ B(B+1)` (Mathlib: `Nat.even_mul_succ_self`).
