# Three-row shapes `(a,b,3)`: `J*` is an XOR-box in `⟨2,4⟩`, so `|J*|∈{1,2,4}` is even (d = 4)

**Date:** 2026-06-14 (prove session)
**Status:** Sub-family `c=3` (three rows, last row three boxes) — **the first family where the
full two-generator box `{0,2,4,6}` is real** (`|J*|=4`). The whole interior *mechanism* is pinned
down rigorously: a **closed form** for `M_j` (sextic numerator), the **Prop-2 Kummer formula**
(identical in shape to `c=1,2`), the **two-generator structural identity**
`Q3 = (a−1)(b−2)·H − 6!·C(j,6)`, and the **offset theorem** (forced descent, proved in full). The
two Compensation inequalities — the `a`-even one and the `a`-odd two-generator one — are reduced to
clean explicit `2`-adic statements and **verified `m ≤ 79`**; the full hand-proof of those two
inequalities, and the boundary residual, are the named gaps (§7).

Companion to `2026-06-13-threerow-c1-Jstar-even.md` and `2026-06-13-threerow-c2-Jstar-even.md`. In
`c=1` the offset was `j₀∈{0,3}` (one generator); in `c=2` both offsets were `0` but the second
generator `4` debuted (`|J*|≤2`). `c=3` is where the two phenomena **combine**: the offset is again
`j₀∈{0,3}` (parity of `a`), and at the offset-`3` branch **both generators fire at once**, giving
the first honest `|J*|=4`.

---

## 0. Setup and result

Fix `m ≥ 1`, `λ ⊢ 2m`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`,

> `M_j = ⟨s_λ , h₁^{2(m−j)} e₂^j⟩ ∈ ℤ_{≥0}`,  `M₀ = f^λ`,
> `val(j) = j + 2 v₂( C(m,j) M_j )` (`= ∞` when `C(m,j)M_j = 0`),
> `V = min_j val(j)`,  `J* = { j : val(j) = V }`.

Leading-π dichotomy (proved earlier, Lean-checked; `π = 1+i`): `C ≡ |J*| (mod π)`, so
**`|J*|` odd ⟹ `G_λ(i) ≠ 0`**; the ties are exactly the shapes with `|J*|` even. Throughout `s₂`
is the binary digit sum, `v₂` the 2-adic valuation, `v₂C(N,k)=s₂(k)+s₂(N−k)−s₂(N)` (Kummer),
`C(n,k)=0` for `k<0` or `k>n`, and `val(j) ≡ j (mod 2)`, so **only `j ≡ j₀ (mod 2)` can tie** with
the minimizer `j₀`.

**Main result.** For every three-row shape `λ = (a,b,3) ⊢ 2m` (`a ≥ b ≥ 3`, `a+b = 2m−3`). Since
`a+b` is odd, exactly one of `a,b` is even, and the parity of `a` fixes the offset `j₀ = min J*`:

> **Theorem.**
> - **`a` even** (so `b` odd): `j₀ = 0` and `J*` is one of `{0}`, `{0,2}`, `{0,4}`:
>   - `J* = {0,4}` ⟺ `a ≡ 0` and `b ≡ 1 (mod 4)`   (generator `4`);
>   - `J* = {0,2}` ⟺ `a ≡ 2` and `b ≡ 3 (mod 4)`   (generator `2`);
>   - `J* = {0}` otherwise.  In all cases `|J*| ≤ 2`.
> - **`a` odd** (so `b` even): `j₀ = 3` and `J*` is one of `{3}`, `{3,5,7,9}`:
>   - `J* = {3,5,7,9} = 3 + ⟨2,4⟩` ⟺ `a ≡ 1` and `b ≡ 2 (mod 4)`   (**both** generators, `|J*|=4`);
>   - `J* = {3}` otherwise.
>
> In every case `J*−j₀` is an XOR-closed `2`-adic box inside `⟨2,4⟩ = {0,2,4,6}`, so
> **`|J*| ∈ {1,2,4}` is even on every tie.** The two generators are *independent* in the `a`-even
> branch (never both: the tie classes `(0,1)` and `(2,3)` are disjoint mod `4`) and **locked
> together** in the `a`-odd branch (both fire ⟺ `a≡1, b≡2`).

The minimal full-box witness is `λ = (9,6,3) = 3·(3,2,1)`, `m = 9`, `J* = {3,5,7,9}` — `a=9≡1`,
`b=6≡2 (mod 4)`. *(Theorem verified against the Murnaghan–Nakayama definition for all interior-tie
`λ=(a,b,3)`, `m ≤ 39`: 630 shapes, 0 mismatch; `|J*|` distribution `{1:405, 2:153, 4:72}`.)*

---

## 1. The closed form for `M_j`

Work in `N=3` alphabet variables `s,t,u`: `h₁=e₁=s+t+u`, `e₂=st+su+tu`. With
`D_j(p,q,r) = [s^p t^q u^r]\,e₁^{2(m−j)}e₂^j` and the three-row Jacobi–Trudi expansion (Prop. 1 of
the `c=1` note), `λ=(a,b,3)` gives (Jacobi–Trudi index `c−3=0`, all six `S₃` terms survive with
`r∈{1,2,3}`):

> `M_j = D_j(a,b,3) − D_j(a,b+1,2) − D_j(a+1,b−1,3) + D_j(a+1,b+1,1) + D_j(a+2,b−1,2) − D_j(a+2,b,1).`

Referencing each `D_j` to the **`t`-coordinate** (the `q`-value, near `b`), the trinomial sums
collapse to a single base binomial times a rational function:

> **Lemma 1 (`c=3` closed form).** For `b ≥ 3` and interior `0 ≤ j ≤ b`:
> `  M_j = \dfrac{C(N,\,b−j)\,(a−b+1)\,Q_3(a,b,j)}{6\,(a+4−j)(b+1−j)(b+2−j)(b+3−j)}`,  `N = 2(m−j)`,
> where `Q_3` is the **sextic-in-`j`**
> `  Q_3(a,b,j) = (a−1)(a+3)(a+4)(b−2)(b+2)(b+3)
>              − 6(a−1)(a+3)(b−2)(b+2)\,C(j,1)
>              − 6(a−1)(b−2)(ab+a+2b)\,C(j,2)
>              + 36(a−1)(b−2)\,C(j,3) + 72(a−1)(b−2)\,C(j,4) − 720\,C(j,6),`
> and `Q_3(a,b,0) = (a−1)(a+3)(a+4)(b−2)(b+2)(b+3) =: K_3`.

The denominator `6·(a+4−j)·∏_{i=1}^{3}(b+i−j)` and the `j`-free prefactor `(a−b+1)` are the exact
`c`-graded analogues of the `c=2` shape `2·(a+3−j)·∏_{i=1}^{2}(b+i−j)`. *(Lemma 1 verified vs
Murnaghan–Nakayama for all `λ=(a,b,3)`, `m ≤ 24`: 2170 cases, 0 mismatch.)*

Because `a+b=2m−3` is **odd**, `a−b+1` is **even**, but it is `j`-free and so cancels in every
`val`-difference, exactly as before.

### The two-generator structural identity

Collecting `Q_3` in the binomial basis `C(j,k)` reveals that **every coefficient except the top one
carries the factor `(a−1)(b−2)`**:

> **Identity.** `  Q_3(a,b,j) = (a−1)(b−2)\,H(a,b,j) \ − \ 6!\,C(j,6),`  where
> `  H = (a+3)(a+4)(b+2)(b+3) − 6(a+3)(b+2)C(j,1) − 6(ab+a+2b)C(j,2) + 36\,C(j,3) + 72\,C(j,4).`

This is the precise `c=3` analogue of the `c=2` decomposition `Q = a(b−1)(W−2j²) + 4!\,C(j,4)`. The
**inhomogeneous tip** is `−6!\,C(j,6) = −P_6(j)`, `P_6(j)=j(j−1)\cdots(j−5)`, vanishing at the
**six** points `j∈\{0,…,5\}` and first nonzero at `P_6(6)=720` — the seed of the *second* generator.
The **heavy factor** `(a−1)(b−2)` is

> **odd** when `a` is even (`a−1,b−2` both odd), and **divisible by `4`** when `a` is odd
> (`a−1,b−2` both even).

This single fact is what splits the theorem into its two branches. *(Identity checked symbolically:
`Q_3 − [(a−1)(b−2)H − 720C(j,6)] ≡ 0`.)*

---

## 2. The Prop-2 Kummer formula

The numerator falling factorials `(b)_j,(a+3)_j` of `C(2m,b)/C(N,b−j) = (2m)_{2j}/[(b)_j(a+3)_j]`
merge with the four denominator linear factors of Lemma 1, via the telescopes
`(b)_j/[(b+1−j)(b+2−j)(b+3−j)] = (b+3)_j/[(b+1)(b+2)(b+3)]` and `(a+3)_j/(a+4−j)=(a+4)_j/(a+4)`,
into `v₂((b+3)_j)` and `v₂((a+4)_j)`. With `v₂((2m)_{2j})=j+v₂((m)_j)` and `C(m,j)=(m)_j/j!` the
`(m)_j` cancels and one obtains the **same shape as `c=2`** with `(a+3,b+2)→(a+4,b+3)`:

> **Proposition 2 (`c=3`).** For `1 ≤ j ≤ b`,
> `  Δ(j) := val(j) − val(0) = j − 2 s₂(j) + 2 v₂C(a+4,\,j) + 2 v₂C(b+3,\,j)
>          + 2\big[\, v₂Q_3(a,b,j) − v₂Q_3(a,b,0) \,\big].`

*(Verified against direct `val(j)−val(0)`, all `λ=(a,b,3)`, `m ≤ 21`: 1275 cases, 0 violations.)*

Write `T_3(j) := v₂C(a+4,j) + v₂C(b+3,j) + v₂Q_3(j) − v₂Q_3(0)`, so `Δ(j) = j − 2s₂(j) + 2T_3(j)`.
The **skeleton** `S(j)=j−2s₂(j)+2v₂C(a+4,j)+2v₂C(b+3,j)` is the proved two-row `Δ` of the shape
`(a+3,b+3)`; the correction `2[v₂Q_3(j)−v₂Q_3(0)]` is where both generators live.

---

## 3. The offset theorem (forced descent) — **proved in full**

> **Proposition 3.** `Δ(1)=1` when `a` is even; `Δ(1)=−1, Δ(2)=−2, Δ(3)=−3` when `a` is odd.
> Hence `j₀ = min J* = 0` for `a` even and `j₀ = 3` for `a` odd.

*Proof.* Use `Q_3(j)=(a−1)(b−2)H(j) − 720C(j,6)` and note `C(j,6)=0` for `j ≤ 5`, so
`Q_3(j)=(a−1)(b−2)H(j)` there; the heavy factor only shifts `v₂` by a `j`-independent constant that
cancels in `v₂Q_3(j)−v₂Q_3(0)`. Thus for `j ≤ 5`, `v₂Q_3(j)−v₂Q_3(0)=v₂H(j)−v₂H(0)`,
`H(0)=(a+3)(a+4)(b+2)(b+3)`.

**`a` even (`b` odd).** Then `a+3,b+2` odd, `a+4,b+3` even; `v₂H(0)=v₂(a+4)+v₂(b+3)`.
`H(1)=(a+3)(b+2)[(a+4)(b+3)−6]`: `(a+3)(b+2)` odd, `(a+4)(b+3)≡0 (mod 4)` so `(a+4)(b+3)−6≡2`,
giving `v₂H(1)=1`. Also `v₂C(a+4,1)=v₂(a+4)`, `v₂C(b+3,1)=v₂(b+3)`. Then
`Δ(1)=1−2+2[v₂(a+4)+v₂(b+3)+v₂H(1)−v₂H(0)]=−1+2·1=1`.

**`a` odd (`b` even).** Then `a+3,b+2` even, `a+4,b+3` odd; `v₂H(0)=v₂(a+3)+v₂(b+2)`.
- `j=1`: `H(1)=(a+3)(b+2)[(a+4)(b+3)−6]` with `(a+4)(b+3)` odd ⟹ `[\cdots]` odd ⟹
  `v₂H(1)=v₂(a+3)+v₂(b+2)=v₂H(0)`; `v₂C(·,1)=0`; so `Δ(1)=1−2+0=−1`.
- `j=2`: `v₂C(a+4,2)=v₂(a+3)−1`, `v₂C(b+3,2)=v₂(b+2)−1` (odd·even/2). `H(2)=(a+3)(b+2)[(a+4)(b+3)
  −12]−6(ab+a+2b)`; the first summand has `v₂≥2`, the second `−6·(\text{odd})` has `v₂=1`
  (`ab+a+2b` is odd for `a` odd, `b` even), so `v₂H(2)=1`. Then
  `Δ(2)=2−2+2[(v₂(a+3)−1)+(v₂(b+2)−1)+1−v₂(a+3)−v₂(b+2)]=−2`.
- `j=3`: `v₂C(a+4,3)=v₂(a+3)−1`, `v₂C(b+3,3)=v₂(b+2)−1`. `H(3)=(a+3)(b+2)[(a+4)(b+3)−18]
  −18(ab+a+2b)+36`: summands have `v₂≥2, =1, =2` respectively, so `v₂H(3)=1`. Then
  `Δ(3)=3−2·2+2[(v₂(a+3)−1)+(v₂(b+2)−1)+1−v₂(a+3)−v₂(b+2)]=−1−2=−3`.

In each `a`-odd line the unique `v₂=1` summand of `H(j)` forces `v₂H(j)=1` (strict min), so the
computation is rigorous. ∎

So for `a` even the analysis is anchored at `j₀=0`; for `a` odd it is a **forced descent**
`val(0) > val(1) > val(2) > val(3)` anchored at `j₀=3`, exactly as in the `c=1` `a`-odd family.

---

## 4. The `a`-even branch: `|J*| ≤ 2` (one generator)

> **Compensation Lemma A (`a` even).** For `a` even, `b` odd, `1 ≤ j ≤ b`:
> `  T_3(j) \ \ge\ 1 − v₂(j),`  hence `Δ(j) ≥ B(j) := j − 2s₂(j) + 2 − 2v₂(j)`.

*Reduction.* Since `a` even ⟹ `a+4` even (`φ₁:=v₂(a+4)`) and `b+3` even (`φ₂:=v₂(b+3)`), with
`v₂Q_3(0)=φ₁+φ₂` (heavy factor odd). From `j\,C(a+4,j)=(a+4)\,C(a+3,j−1)` (`a+3` odd):
`v₂C(a+4,j)+v₂(j)=φ₁+v₂C(a+3,j−1) ≥ φ₁`. Therefore
`T_3(j)+v₂(j) ≥ φ₁ + v₂C(b+3,j) + v₂Q_3(j) − φ₁ − φ₂ = v₂C(b+3,j) + v₂Q_3(j) − φ₂`, and the lemma
reduces to

> **(L3″-A)** `  v₂C(b+3,\,j) + v₂Q_3(j) \ \ge\ v₂(b+3) + 1 \qquad (a \text{ even}).`

*(Both this and the `(a+4)`-symmetric variant verified: all `a`-even `λ`, `1≤j≤b`, `m ≤ 79`; 0
failures, tight `12114` times.)* **(L3″-A) is the one unproved interior inequality of the `a`-even
branch — see §7.** The naive route `v₂Q_3(j) ≥ v₂(j)+1` is **false** at `j=4` (`v₂(j)=2`): there the
binomial slack `v₂C(a+3,j−1)` must pay the deficit. This is precisely the `c=2` Number-Lemma
phenomenon, now carried by the tip `−6!C(j,6)`; the missing step is the subset identity
`C(F+s,j)C(j,6)=C(F+s,6)C(F+s−6,j−6)` analysis with the *odd* heavy factor.

*From the lemma to the box.* `B(j)=j+2−2s₂(j)−2v₂(j)`. By the `c=2` Lemma B (proved there),
`B(j)≥0` with `B(j)=0` **iff `j∈{2,4}`** and `B(j)≥1` otherwise. Hence `Δ(j)≥0` for all interior
`j≥1`, with equality only possible at `j∈{2,4}`; combined with the boundary bound (§6) this gives
`0∈J*` and `J*⊆{0,2,4}`. The two ties are governed by the heavy-factor-free valuations
`v₂H(2),v₂H(4)`:
- `Δ(2)=2[v₂H(2)−2]`, so `Δ(2)=0 ⟺ v₂H(2)=2 ⟺ a≡2, b≡3 (mod 4)`;
- `Δ(4)=0 ⟺ a≡0, b≡1 (mod 4)`.

These mod-`4` classes are **disjoint**, so `j=2` and `j=4` never tie simultaneously: `|J*| ≤ 2`.
*(The two tie-congruences verified, all `a`-even `λ`, `m ≤ 39`.)*

---

## 5. The `a`-odd branch: `|J*| ∈ {1,4}` (both generators, locked)

Here the heavy factor `(a−1)(b−2)` is divisible by `4`, and `a+4,b+3` are **odd**. After the forced
descent (Prop. 3) the baseline is `j₀=3`; set `\tildeΔ(j) := val(j) − val(3) = Δ(j) + 3`.

> **Compensation Lemma B (`a` odd).** For `a` odd, `b` even, `4 ≤ j ≤ b`:
> `  \tildeΔ(j) = j + 3 − 2 s₂(j) + 2\,U(j) \ \ge\ 0, \qquad
>     U(j) := v₂C(a+4,j) + v₂C(b+3,j) + v₂Q_3(j) − v₂Q_3(0),`
> with equality (a tie at `j`) possible only at `j ∈ {5,7,9}`, and then **all three tie together**:
> `\tildeΔ(5)=\tildeΔ(7)=\tildeΔ(9)=0 ⟺ a ≡ 1` and `b ≡ 2 (mod 4)`.

*(Verified: `\tildeΔ(j) ≥ 0` for all `a`-odd `λ`, `4 ≤ j ≤ b`, `m ≤ 79`; 0 failures. The tie pattern
— all-or-nothing on `{5,7,9}`, congruence `a≡1,b≡2` — verified `m ≤ 39`.)*

This is the genuinely new content: the **box is the full group `⟨2,4⟩`**. `\tildeΔ(5)=0` is the
generator-`2` tie (relative to `3`), `\tildeΔ(9)=0` (i.e. `9 = 3+6`, `6 = 2\oplus4`) is forced by the
XOR-closure, and `\tildeΔ(7)=0` is generator-`4`. The **lock** (`|J*|∈{1,4}`, never `2`) is the
statement that in the `a`-odd branch the two generators share the single switch `a≡1, b≡2 (mod 4)` —
a sharper coupling than the `a`-even branch, where they are independent. The full hand-proof of
Compensation Lemma B (the floor is *not* a single clean `c−v₂(j)`; the per-`v₂(j)` minima of `U` are
`{0:−6,1:−5,2:−5,3:−4,4:−5,…}`) is the **two-generator Number Lemma** still open — see §7.

---

## 6. Boundary regime and the residual

A tie point is *interior* iff it is `≤ b`. For `a` even the tie points `{2,4}` are interior once
`b ≥ 4` (and `b≥3` carries only `j=2`); for `a` odd the box reaches `j₀+6 = 9`, interior once
`b ≥ 9`. For smaller `b` the upper tie carriers fall on the boundary `j∈\{b+1,b+2,b+3\}`, where (by
Lemma 1) only `M_{b+1},M_{b+2},M_{b+3}` are nonzero and exactly the ones of the right parity can tie.

As in `c=1,2`, the comparison `val(b+i) > val(j₀)` (for the non-tie boundary points) crosses from the
factored interior form to the boundary values and we do **not** have a hand proof. It is verified
**`m ≤ 79`** together with the small-`b` boundary-tie families (which match the Theorem's
congruences). This is the standard lone residual of the family.

---

## 7. What is proved, what is verified, and the named gaps

**Proved in full (symbolic/closed-form, no `m`-cap):**
- Lemma 1 closed form (numerator `Q_3`, denominator, prefactor) — and its verification harness.
- The two-generator structural identity `Q_3=(a−1)(b−2)H − 6!C(j,6)`.
- Proposition 2 (`Δ` Kummer formula) shape.
- **Proposition 3 (offset theorem / forced descent)** — `j₀∈\{0,3\}` by parity of `a`, fully proved.
- The **reductions**: Compensation Lemma A ⟸ (L3″-A); the box `J*⊆\{0,2,4\}` and `|J*|≤2` for `a`
  even ⟸ Compensation Lemma A + the disjoint mod-`4` tie classes; the `a`-odd box `⟨2,4⟩` ⟸
  Compensation Lemma B.

**Verified computationally (not yet hand-proved), with `m`-range stated:**
- Lemma 1 vs Murnaghan–Nakayama (`m≤24`); Prop. 2 (`m≤21`); full Theorem (`m≤39`).
- **(Gap 1) (L3″-A):** `v₂C(b+3,j)+v₂Q_3(j) ≥ v₂(b+3)+1` for `a` even (`m≤79`). Closing it needs the
  `C(j,6)`-tip Number Lemma with an *odd* heavy factor.
- **(Gap 2) Compensation Lemma B:** `\tildeΔ(j)≥0` for `a` odd, `4≤j≤b` (`m≤79`), and the
  all-or-nothing tie pattern on `{5,7,9}`. This is the **two-generator Number Lemma** — the real
  prize: it must produce *two* simultaneous deficits (`v₂(j)=1` and `v₂(j)=2`) keyed on the single
  switch `a≡1,b≡2 (mod 4)`.
- **(Gap 3) Boundary residual** `val(b+i) > val(j₀)` (`m≤79`), the standard lone residual.

**Net:** the `c=3` closed form and the two-generator decomposition are derived and verified; the
offset structure is proved; even-`|J*|` is **reduced to two explicit `2`-adic inequalities** (Gaps
1–2), each verified `m≤79`. This is the on-ramp to the general `e₂ mod 2` wall: the engine
`Q_{(a,b,c)} = (\text{heavy})·H − (2c)!\,C(j,2c)` and the offset-by-`parity(a)` law now both visibly
`c`-graded.

---

## 8. Diagnosis — how `|J*|=4` is born

`c=1`: tip `−j(j−1)` (root `\{0,1\}`), one generator `2`, offset `j₀∈\{0,3\}` by parity of `a`.
`c=2`: tip `+4!\,C(j,4)` (root `\{0,…,3\}`), generators `2` **or** `4` but never both (`|J*|≤2`),
offset `0`. `c=3` **superposes** the two: the tip is `−6!\,C(j,6)` (root `\{0,…,5\}`, first nonzero
`720` at `j=6`), and the offset is again `j₀∈\{0,3\}`. The new phenomenon is that **at the
offset-`3` (`a`-odd) branch the heavy factor `(a−1)(b−2)` is divisible by `4`**, which is exactly the
arithmetic that lets *both* the `v₂(j)=1` deficit (generator `2`, `j=5`) and the `v₂(j)=2` deficit
(generator `4`, `j=7`) be paid at once — and their XOR `j=9` for free. So the full box `⟨2,4⟩`
(`|J*|=4`) is the `a`-odd branch's two-generator coincidence, switched on by `a≡1, b≡2 (mod 4)`. The
general pattern is now legible: `c` controls the tip order `(2c)!\,C(j,2c)` (generators
`2,4,…,2c`), and `parity(a)` controls the offset; coexistence of generators is a heavy-factor
divisibility event at the shifted offset.

### Files
`projects/code/threerow-c3/`: `c3factor.py` (Lemma 1 derivation), `c3verify.py` (Lemma 1 vs MN +
binomial-basis check), `c3prop2.py` (Prop. 2), `c3decomp.py` (structural identity),
`c3comp.py`,`c3reduc.py`,`c3lemtest.py` (Compensation Lemmas A/B, reductions, `m≤79`),
`c3census.py`,`c3detail.py` (box by class), `c3full.py` (full Theorem, `m≤39`). Shared: `mn.py`.
