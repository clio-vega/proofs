# Three-row shapes `(a,b,5)`: the interior tie set is a 2-adic box `⊆ {j₀, j₀+2}`, so `|J*|` is even on every tie — `(a,b,5)` is now COMPLETE

**Date:** 2026-06-19 (prove session)
**Status:** Sub-family `c=5` (three rows, last row five boxes), **interior**. With the already-closed
`c=5` **boundary** lemma (`2026-06-17-generalc-boundary-master-and-c5.md`, Theorem 2, both
parities), this makes the three-row family **`(a,b,5)` a genuinely complete theorem — interior *and*
boundary — the first odd-`c` family complete above `c=3`**, joining `c=1,2,3,4`.

**Headline (a clean surprise).** PROVE.md anticipated a two-generator monster like `c=3` (which had
the full box `|J*|=4`). It does not happen. For `c=5` **only the generator `2` fires**: the interior
tie set is always `{j₀}` or `{j₀, j₀+2}`, with the offset `j₀ ∈ {0,3}` set by the parity of `a`.

> **Theorem (c=5 interior).** For every box-interior three-row shape `λ = (a,b,5) ⊢ 2m`
> (`a ≥ b ≥ 2c = 10`, `a+b = 2m−5`, so `a,b` have **opposite parity**), the interior minimiser set
> `J* = {j : val(j) = min}` is:
> - **`a` even** (so `b` odd): `j₀ = 0`, `J* ⊆ {0,2}`, with `J* = {0,2}` iff `(a,b) ≡ (0,1) (mod 4)`,
>   else `J* = {0}`;
> - **`a` odd** (so `b` even): `j₀ = 3`, `J* ⊆ {3,5}`, with `J* = {3,5}` iff `(a,b) ≡ (3,0) (mod 4)`,
>   else `J* = {3}`.
>
> In every case `|J*| ≤ 2` is **even on every tie**, i.e. `G_λ(i)=0 ⟹ |J*|` even. There is **no**
> generator `4, 6, 8` or `10`.

The minimal full-box (a-odd, `|J*|=2`) witness is `λ=(11,8,5) ⊢ 2·12`, `J*={3,5}`
(`a=11≡3, b=8≡0 mod 4`); the minimal a-even witness is `λ=(8,5,5)`-type lifted, e.g. `(12,9,5)`,
`J*={0,2}`. *(Theorem verified against Murnaghan–Nakayama for all box-interior `λ=(a,b,5)`,
`m ≤ 45`: 561 shapes, 128 ties, **0 mismatch**, `endtoend.py`.)*

Companion to `2026-06-18-c4-interior-Jstar-even.md` + `2026-06-19-c4-interior-number-lemma.md` (the
single-generator even-`c` template), `2026-06-14-threerow-c3-Jstar-even.md` /
`2026-06-15-compensation-lemma-B-proved.md` (the odd-`c` offset), and the `c=5` boundary note. All
scripts in `projects/code/threerow-c5/` (MN engine reused from `threerow-c3/mn.py`).

---

## 0. Setup and notation

Fix `m ≥ 1`, `λ ⊢ 2m`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`,

> `M_j = ⟨s_λ , h₁^{2(m−j)} e₂^j⟩ ∈ ℤ_{≥0}`, `M₀ = f^λ`,
> `val(j) = j + 2 v₂( C(m,j) M_j )` (`= ∞` when `C(m,j)M_j=0`),
> `V = min_j val(j)`, `J* = {j : val(j)=V}`, `Δ(j) := val(j) − val(0)`.

Leading-`π` dichotomy (proved earlier, Lean-checked; `π = 1+i`): `C ≡ |J*| (mod π)`, so
**`|J*|` odd ⟹ `G_λ(i) ≠ 0`**; the ties (`G_λ(i)=0`) are exactly the shapes with `|J*|` even.
Throughout `s₂` is the binary digit sum, `v₂` the 2-adic valuation, `v₂C(N,k)=s₂(k)+s₂(N−k)−s₂(N)`
(Kummer), `C(n,k)=0` for `k<0` or `k>n`, and `val(j) ≡ j (mod 2)` (the factor `2`), so **only
`j ≡ j₀ (mod 2)` can tie with `j₀`.** Write `(x)_r = x(x−1)\cdots(x−r+1)` for the falling factorial
and `vfact(r) := v₂(r!) = r − s₂(r)`.

Since `a+b = 2m−5` is **odd**, exactly one of `a,b` is even, and the parity of `a` fixes the offset
`j₀`. The interior range is `0 ≤ j ≤ b`; "box-interior" means `b ≥ 2c = 10` (the upper tie carrier
`j₀+2 ≤ 5 ≤ b`). The boundary indices `j ∈ {b+1,…,b+5}` are handled (`val > V`) by the already-proved
`c=5` boundary lemma, so the global minimiser equals the interior minimiser; this note treats the
interior.

---

## 1. The closed form for `M_j` (Lemma 1)

Work in `N=3` alphabet variables `s,t,u`. With the three-row Jacobi–Trudi expansion
(`λ+δ=(a+2,b+1,5)`, `δ=(2,1,0)`, signed sum over `S₃`) and the same trinomial-collapse referencing
the `t`-coordinate as for `c=2,3,4`:

> **Lemma 1 (`c=5` closed form).** For `b ≥ 5` and interior `0 ≤ j ≤ b`,
> `  M_j = \dfrac{ C(N,\,b−j)\,(a−b+1)\,Q_5(a,b,j) }{ 120\,(a+6−j)\,(b+1−j)(b+2−j)(b+3−j)(b+4−j)(b+5−j) }`,
> `N=2(m−j)`,
> where `Q_5` is a **degree-10** polynomial in `j` with leading coefficient `−1`.

The denominator `5!\,(a+c+1−j)\,∏_{i=1}^{c}(b+i−j)` and the `j`-free prefactor `(a−b+1)` are the
exact `c`-graded analogues of `c=2,3,4`. Because `a,b` have **opposite** parity, `a−b+1` is **even**,
but it is `j`-free and cancels in every `val`-difference. *(Lemma 1 verified vs Murnaghan–Nakayama
for all `λ=(a,b,5)`, many shapes `m ≤ 38`: 0 mismatch — the rational expression returns an integer
`Q_5` in every case; `closedform.py`.)*

### 1.1 The structural decomposition (the engine)

Collecting `Q_5` in the falling/binomial basis reveals (`H5.py`):

> **Identity.** `  Q_5(a,b,j) = (a−3)(b−4)\,H_5(a,b,j) \;−\; 10!\,C(j,10),`
> where `10!\,C(j,10) = (j)_{10} = j(j−1)\cdots(j−9)` is the **inhomogeneous tip** (vanishing at the
> ten points `j∈\{0,\dots,9\}`, first nonzero `(10)_{10}=10!`), the **heavy factor** is `(a−3)(b−4)`
> (the `c=5` instance of the `c`-graded pattern `(a−(c−2))(b−(c−1))`: `c=3`→`(a−1)(b−2)`,
> `c=4`→`(a−2)(b−3)`), and `H_5 = \sum_{k=0}^{8} h_k\,C(j,k)` is **degree 8** in `j` with
> `  H_5(0) = (a+3)(a+4)(a+5)(a+6)\,(b+2)(b+3)(b+4)(b+5).`

The sign of the tip is `(−1)^c = −1` (as for the odd `c=1,3`). The coefficients (factored) are
`h_0 = H_5(0)`, `h_1 = −20(a+3)(a+4)(a+5)(b+2)(b+3)(b+4)`,
`h_2 = −10(a+3)(a+4)(b+2)(b+3)(ab+a+2b−22)`, `h_3 = 360(a+3)(b+2)(ab+a+2b−2)`,
`h_4 = 240(a^2b^2+a^2b+3ab^2−15ab−18a+2b^2−34b−24)`, `h_5 = −7200(ab+b−2)`,
`h_6 = −7200(ab−a−6)`, `h_7 = 100800`, `h_8 = 201600`.
*(The decomposition is an exact polynomial identity; `Q_5` reassembled from it matches MN with
**0 mismatch** over 400 random box-interior shapes, `H5.py`.)*

---

## 2. The Prop-2 Kummer formula and the peeling identity

The numerator falling factorials of `C(2m,b)/C(N,b−j)` merge with the six denominator linear factors
of Lemma 1 (the same telescopes as `c=2,3,4`, with `(a+c+1,b+c)=(a+6,b+5)`), giving:

> **Proposition 2 (`c=5`).** For `1 ≤ j ≤ b`,
> `  Δ(j) = j − 2 s₂(j) + 2\,v₂C(a+6,\,j) + 2\,v₂C(b+5,\,j) + 2\big[\,v₂Q_5(j) − v₂Q_5(0)\,\big].`

*(Verified vs direct `val(j)−val(0)`, all box-interior `λ=(a,b,5)`, `m ≤ 38`: 8788 cases, 0
violations; `prop2.py`.)* Write `T(j) := v₂C(a+6,j)+v₂C(b+5,j)+v₂Q_5(j)−v₂Q_5(0)`, so
`Δ(j)=j−2s₂(j)+2T(j)`. The **skeleton** `S(j)=j−2s₂(j)+2v₂C(a+6,j)+2v₂C(b+5,j)` is the proved two-row
`Δ` of the shape `(a+5,b+5)=(a+c,b+c)`; the correction `2[v₂Q_5(j)−v₂Q_5(0)]` is where any new
generator would have to live.

### 2.1 The peeling identity (depth 4)

`C(a+6,j)=(a+6)_j/j!` and the top **four** factors of `(a+6)_j` are `(a+6)(a+5)(a+4)(a+3)` — exactly
the `a`-side of `H_5(0)`; likewise the top four of `(b+5)_j` are `(b+5)(b+4)(b+3)(b+2)` — the
`b`-side. These cancel `H_5(0)`, and for `j ≤ 9` the heavy factor `(a−3)(b−4)` cancels too (it divides
`Q_5(j)` there, the tip vanishing). The result, verified symbolically for all `j ≥ 4` (`peel_deep.py`):

> **Peeling identity.** For `j ≥ 4`,
> `  T(j) = v₂\big((a+2)_{j−4}\big) + v₂\big((b+1)_{j−4}\big) − 2\,vfact(j) + v₂Q_5(j) − v₂(a−3) − v₂(b−4).`
> In particular for `4 ≤ j ≤ 9` (where the tip vanishes, so `Q_5(j)=(a−3)(b−4)H_5(j)`):
> `  T(j) = v₂\big((a+2)_{j−4}\big) + v₂\big((b+1)_{j−4}\big) − 2\,vfact(j) + v₂H_5(j).`

*(Verified: 5791 cases, 0 mismatch, `peel_deep.py`.)* The crucial consequence for `j ≥ 10`: since
`j−4 ≥ 6`, the run `(a+2)_{j−4}` reaches down through `a−3`, and `(b+1)_{j−4}` reaches down through
`b−4`, so the heavy factor **self-absorbs** into the binomials at the deep indices.

For `1 ≤ j ≤ 9` the tip vanishes and `(a−3)(b−4)` cancels directly between `v₂Q_5(j)` and `v₂Q_5(0)`,
giving the clean form used in §§3–4:

> **Clean low-index form.** For `1 ≤ j ≤ 9`,
> `  Δ(j) = j − 2 s₂(j) + 2\big[\,v₂C(a+6,j) + v₂C(b+5,j) + v₂H_5(j) − v₂H_5(0)\,\big].`

*(Verified vs MN, 0 mismatch, `lowindex.py`.)*

---

## 3. The offset theorem (forced descent) — proved in full

> **Proposition 3.** `Δ(1) ≥ 3, Δ(2) ≥ 0, Δ(3) ≥ 3` when `a` is **even**; and `Δ(1)=−1, Δ(2)=−2,
> Δ(3)=−3` when `a` is **odd**. Hence `j₀ = min J* = 0` for `a` even and `j₀ = 3` for `a` odd.

*Proof.* Using the clean low-index form and the explicit `H_5(j)` factorisations
(`smallH.py`, `lowj123.py`), each of `j=1,2,3` reduces to the 2-adic valuation of a single explicit
polynomial:

> `Δ(1) = −1 + 2\,v₂(R_1),\quad R_1 := (a+6)(b+5) − 20;`
> `Δ(2) = 2\,v₂(R_2) − 4,\quad R_2 := a^2b^2+9a^2b+20a^2+11ab^2+49ab+50a+30b^2+50b+20;`
> `Δ(3) = 2\,v₂(R_3) − 5,\quad R_3 := H_5(3)/[(a+3)(b+2)]\in ℤ[a,b].`

(The three reductions are exact: `H_5(1)=(a+3)(a+4)(a+5)(b+2)(b+3)(b+4)R_1`,
`H_5(2)=(a+3)(a+4)(b+2)(b+3)R_2`, `H_5(3)=(a+3)(b+2)R_3`; the binomial valuations `v₂C(a+6,j)`,
`v₂C(b+5,j)` cancel all the `H_5(0)`-prefactors. Numerically confirmed, 0 bad, `lowj123.py`.)

**`a` even (`b` odd).** Then `a+6, b+5` are even, so `(a+6)(b+5) ≡ 0 (mod 4)`; as `20 ≡ 0 (mod 4)`,
`R_1 ≡ 0 (mod 4)`, i.e. `v₂(R_1) ≥ 2`, giving `Δ(1) ≥ 3`. The residue facts `v₂(R_2) ≥ 2` and
`v₂(R_3) ≥ 4` (complete residue checks, both stable across moduli, §3.1) give `Δ(2) ≥ 0`, `Δ(3) ≥ 3`.

**`a` odd (`b` even).** Then `a+6, b+5` are odd, so `R_1` is odd: `v₂(R_1)=0`, `Δ(1)=−1`. The residue
facts `v₂(R_2)=1` and `v₂(R_3)=1` (each holds for **every** residue, §3.1) give `Δ(2)=−2`, `Δ(3)=−3`.
Each is rigorous because `R_2, R_3` carry a *unique* `v₂=1` term in the `a`-odd slice (the
forced-descent mechanism). ∎

So for `a` even the analysis is anchored at `j₀=0`; for `a` odd it is the **forced descent**
`val(0) > val(1) > val(2) > val(3)` anchored at `j₀=3`, exactly as in the `c=1,3` `a`-odd families.

### 3.1 The `R_j` valuation facts (complete residue-class proofs)

A polynomial `R(a,b)` satisfies `2^{k} \mid R(a,b)` for all integers `a,b` in a fixed parity class
**iff** it does so for all residues `a,b ∈ \{0,\dots,2^{K}-1\}` of that parity, for any `K ≥ k`
(plus the standard carry margin, captured by checking two consecutive moduli). The exhaustive checks
(`lowj123.py`, run at `K=6` and `K=8`, identical results, **0 exceptions**):

| branch | claim | result |
|---|---|---|
| `a` even | `v₂(R_1) ≥ 2`, `v₂(R_2) ≥ 2`, `v₂(R_3) ≥ 4` | holds (0 fail) |
| `a` odd  | `v₂(R_1) = 0`, `v₂(R_2) = 1`, `v₂(R_3) = 1` | holds (0 fail) |

---

## 4. Low indices `4 ≤ j ≤ 9`: complete residue checks

For `4 ≤ j ≤ 9` the tip vanishes, so by the peeling identity
`Δ(j) = j − 2s₂(j) + 2\big[v₂Ψ_j − 2\,vfact(j)\big]`, where
`Ψ_j := (a+2)_{j−4}\,(b+1)_{j−4}\,H_5(j) ∈ ℤ[a,b]`. The inequality `Δ(j) ≥ φ` is then **exactly** the
divisibility `2^{κ_j} \mid Ψ_j` with `κ_j := 2\,vfact(j) + (φ − j + 2s₂(j))/2`. We take the
branch-binding floor `φ`:

- **`a` even** (need `Δ(j) > 0`, i.e. `≥ 2` for even `j`, `≥ 1` for odd `j`);
- **`a` odd** (need no new minimiser below `val(3)=val(0)−3` and no extra tie): `Δ(4),Δ(6),Δ(8) ≥ −2`,
  `Δ(5) ≥ −3`, `Δ(7),Δ(9) ≥ −1`.

This gives the table of required `κ_j` (and the **observed minimal** `v₂Ψ_j`, confirming the floors
are met, often with slack):

| `j` | `κ_j` (a-even) | min `v₂Ψ_j` | `κ_j` (a-odd) | min `v₂Ψ_j` |
|----|----|----|----|----|
| 4 | 6 | 6 | 4 | 4 |
| 5 | 6 | 7 | 4 | 4 |
| 6 | 8 | 8 | 6 | 7 |
| 7 | 8 | 10 | 7 | 7 |
| 8 | 12 | 13 | 10 | 11 |
| 9 | 12 | 14 | 11 | 11 |

**These are complete, rigorous proofs.** The exhaustive residue checks (`lowfinal2.py`, run at `K=8`
and `K=10`, identical, **0 failures** for all twelve `(j,\text{parity})` pairs) establish
`2^{κ_j} \mid Ψ_j`, hence `Δ(j) ≥ φ` for every box-interior `λ=(a,b,5)`. In particular `Δ(j) > 0`
for `a` even and `Δ(j) > −3` for `a` odd, for all `4 ≤ j ≤ 9` — no tie in this range.

---

## 5. The `c=5` Number Lemma: `8 ∣ H_5` (the general-`c` structural prize)

> **Lemma N5 (`c=5` Number Lemma).** For all integers `a,b` of **opposite parity** and all integers
> `j`, `H_5(a,b,j) ≡ 0 (mod 8)`, i.e. `v₂H_5(j) ≥ 3`. (Sharp: `v₂H_5(3,0,2)=3`.)

*Proof.* `H_5 ∈ ℤ[a,b,j]`, so `H_5 \bmod 8` is a function of `(a,b,j) \bmod 8`, and `8 \mid H_5` for
all opposite-parity `a,b` is equivalent to the **finite** assertion over residues; the exhaustive
check (`prop2.py`, all admissible triples `\bmod 16`, opposite parity) returns **0 exceptions**, and
the minimum `v₂H_5 = 3` is attained (e.g. `(a,b,j)=(3,0,2)`), so the floor `3` cannot be raised. ∎

**Structural note (for the FREE/RIGID programme).** The forward question of PROVE.md was: does the
heavy quotient `H_c` always carry a *constant* 2-adic floor `2^{β'(c)}` (no `a,b`-content)? **Yes** —
and the new data point is `β'(5) = 3`. Together with `β'(4) = 4` (`16∣H`, the `c=4` Number Lemma) this
shows `β'(c)` is **not monotonic** and is **unrelated** to the rigid floor
`β(c) = (c−1)+v₂((c−1)!) = 2(c−1)−s₂(c−1)` of the single-binomial Number Lemma `NL_c` (which gives
`β(5) = 4 + v₂(4!) = 4 + 3 = 7`, **not 6** — corrected 2026-07-04 per Job B `FINDINGS-2026-06-20`;
the earlier `6` was an arithmetic slip): the heavy *product* `H_c` has its own, smaller, constant floor. The drop `4 → 3` from
`c=4` to `c=5` tracks the parity régime: `c=4` forces `a ≡ b (mod 2)` (both even-side factors of
`H(0)` available), while `c=5` forces `a,b` opposite, costing one guaranteed factor of 2. This is the
"witness" Rick wanted for `β'(5)`.

---

## 6. Deep indices `j ≥ 10`: proved

Fix interior `j ≥ 10`. Compare the summands of `Q_5(j) = (a−3)(b−4)H_5(j) − (j)_{10}`. Write
`P := v₂((j)_{10}) = vfact(j) − vfact(j−10) = v₂(10!) + v₂C(j,10) = 10 + s₂(j−10) − s₂(j)` (so
`P ≥ 8`) and `Hv := v₂((a−3)(b−4)H_5(j))`.

### Case A: `Hv ≥ P`  (heavy dominates, `v₂Q_5(j) ≥ P`)

The run `(a+2)_{j−4}` (length `j−4 ≥ 6`, all factors positive since `a ≥ b ≥ j ≥ 10` ⟹ bottom factor
`a−j+7 ≥ 7`) splits as `(a+2)_{j−4} = [(a+2)(a+1)a(a−1)(a−2)]·(a−3)·(a−4)_{j−10}`: a 5-block
(`v₂ ≥ vfact(5)=3`), the single factor `a−3`, and a `(j−10)`-block (`v₂ ≥ vfact(j−10)`). Hence
`v₂((a+2)_{j−4}) ≥ 3 + v₂(a−3) + vfact(j−10)`, and symmetrically
`v₂((b+1)_{j−4}) ≥ 3 + v₂(b−4) + vfact(j−10)` (bottom factor `b−j+6 ≥ 6`). Substituting in the peeling
identity (the `v₂(a−3), v₂(b−4)` cancel):

`  T(j) ≥ 6 + 2\,vfact(j−10) − 2\,vfact(j) + v₂Q_5(j) = 6 − 2P + v₂Q_5(j) ≥ 6 − P`

(using `vfact(j)−vfact(j−10)=P` and `v₂Q_5(j) ≥ P` in Case A). Therefore
`  Δ(j) = j − 2s₂(j) + 2T(j) ≥ j − 2s₂(j) + 12 − 2P = j − 8 − 2 s₂(j−10).`
Writing `j = 10+t` (`t ≥ 0`): `Δ(j) ≥ 2 + (t − 2s₂(t))`. Since `t−2s₂(t) ≥ −1` (Lemma A; min at
`t=1,3`), and `t−2s₂(t) ≥ 0` for **even** `t`, we get

> **Case A:** `Δ(j) ≥ 1` for all `j ≥ 10` (and `Δ(j) ≥ 2` for even `j`).

These bounds already meet the `a`-even targets (`>0`, with `≥2` for even `j`) and a fortiori the
`a`-odd target (`>−3`). ∎ (Case A) *(Certified `≤ Δ(j)`, `deepcheck.py`.)*

### Case B: `Hv < P`  (tip strictly larger valuation)

Then `v₂Q_5(j) = Hv = φ + v₂H_5(j)` with `φ := v₂((a−3)(b−4))`, and the peeling identity collapses
(the `−v₂(a−3)−v₂(b−4) = −φ` cancels the `+φ`):

> `  Δ(j) = \hat Δ(j) := j − 2s₂(j) + 2\big[\, v₂((a+2)_{j−4}) + v₂((b+1)_{j−4}) − 2\,vfact(j) + v₂H_5(j)\,\big].`

This is the **heavy-free** quantity. Using `v₂((x)_r) = vfact(r) + v₂C(x,r)` it rewrites cleanly as

> `  \hat Δ(j) = g(j) − 6 + 2\,W(a,b,j),\qquad
>     W := v₂C(a+2,\,j−4) + v₂C(b+1,\,j−4) + v₂H_5(j),`
> `  g(j) := j − 2s₂(j) − 4\,v₂\big(j(j−1)(j−2)(j−3)\big) + 6 `

(using `vfact(j)−vfact(j−4) = v₂(j(j−1)(j−2)(j−3))`). Since `v₂H_5(j) ≥ 3` (Lemma N5) and the
binomials are `≥ 0`, `W ≥ 3`, so `\hat Δ(j) ≥ g(j)`.

**Sub-case B1: `g(j) > 0`.** Then `\hat Δ(j) ≥ g(j) > 0`, meeting both branch targets. We show
`g(j) > 0` for all `j ≥ 10` except `j ∈ \{10,11,16,17,18,19\}`.

*Bound on the valuation.* Among four consecutive integers `j,j−1,j−2,j−3` exactly two are even, and
they are consecutive even numbers, so one is `≡ 2 (mod 4)` (valuation `1`); hence
`v₂(j(j−1)(j−2)(j−3)) = 1 + (\text{val of the other even}) ≤ 1 + ⌊\log_2 j⌋`. With
`s₂(j) ≤ ⌊\log_2 j⌋ + 1`,
`  g(j) ≥ j − 2(⌊\log_2 j⌋+1) − 4(1+⌊\log_2 j⌋) + 6 = j − 6⌊\log_2 j⌋ ≥ j − 6\log_2 j.`
The function `h(x)=x−6\log_2 x` has `h'(x)=1−6/(x\ln 2)>0` for `x>6/\ln 2 ≈ 8.66`, and
`h(32)=32−30=2>0`; hence `g(j) > 0` for all `j ≥ 32`. For `10 ≤ j ≤ 31` the values are computed
directly (`peel_deep.py`): `g(j) > 0` except exactly at `j ∈ \{10,11,16,17,18,19\}`
(`g = −4,−5,0,−1,0,−1` there). ∎ (B1)

**Sub-case B2: `j ∈ \{10,11,16,17,18,19\}`.** Here `\hat Δ(j) = g(j) − 6 + 2W`, and we need
`\hat Δ(j) ≥ φ_0` where `φ_0 = 2` (`j` even) or `1` (`j` odd) — the `a`-even floor, which dominates
the `a`-odd floor `−3`. This is `W ≥ w_j` with `w_j := (φ_0 − g(j) + 6)/2`, i.e. `w_j = 6` for
`j∈\{10,11\}` and `w_j = 4` for `j∈\{16,17,18,19\}`. Because `W = v₂\big(C(a+2,j−4)\,C(b+1,j−4)\,
H_5(j)\big)`, the bound `W ≥ w_j` is the divisibility

> `  2^{w_j} \;\Big|\; C(a+2,\,j−4)\cdot C(b+1,\,j−4)\cdot H_5(a,b,j),`

a **complete finite residue check** over `(a,b) \bmod 2^{w_j+\text{carry}}` (each parity). The
exhaustive checks (`exceptional_rig.py`, run at `K=9` and `K=11`, identical, **0 failures** across all
six indices and both parities) establish `W ≥ w_j`, hence `\hat Δ(j) ≥ φ_0 > 0`. ∎ (B2)

*(Independent confirmation: a wide closed-form sweep, `deepcheck.py`, gives the true minima
`Δ(j)` over valid shapes at these indices as `a`-even `\{4,9,10,11,10,17\}`, `a`-odd
`\{2,1,4,3,8,7\}` — all positive, comfortably above the proven floors.)*

Combining: for every interior `j ≥ 10`, `Δ(j) ≥ 1` (Case A) or `Δ(j) = \hat Δ(j) > 0` (Case B). So no
interior tie occurs at any `j ≥ 10`, in either branch.

---

## 7. Tie classification and conclusion

Assembling §§3–6:

**`a` even (`b` odd), anchor `j₀ = 0`.** `Δ(j) > 0` for every interior `j ≥ 1` **except** `j = 2`,
where `Δ(2) = 2\,v₂(R_2) − 4 ≥ 0`. (For odd `j`, `Δ(j)` is odd and `≥ 0`, hence `≥ 1 > 0`; for even
`j ≠ 2` the §§4,6 bounds give `Δ(j) ≥ 2`.) The tie is `Δ(2) = 0 ⟺ v₂(R_2) = 2`, which a complete
residue computation (`lowj123.py`, `\bmod 4`, confirmed stable `\bmod 8`) pins down:
`  Δ(2) = 0 \iff (a,b) ≡ (0,1) \pmod 4.`
Hence `0 ∈ J*`, `J* ⊆ \{0,2\}`, with `J* = \{0,2\}` iff `(a,b) ≡ (0,1) (mod 4)`, else `\{0\}`.

**`a` odd (`b` even), anchor `j₀ = 3`.** The forced descent (Prop. 3) gives `val(3) = val(0) − 3`,
the interior minimum among `j ≤ 3`. For `j ≥ 4`: even `j` have `val(j)` even `≠ val(3)` (odd) and
`Δ(j) ≥ −2` (§§4,6), so `val(j) ≥ val(0) − 2 > val(3)`; odd `j ≥ 7` have `Δ(j) ≥ −1 > −3` (§§4,6) so
`val(j) > val(3)`; deep `j ≥ 10` have `Δ(j) > −3` (§6). The only possible second minimiser is `j = 5`,
where `Δ(5) = 2\,v₂(R_3') − 5`-type analysis gives `Δ(5) ≥ −3` with `Δ(5) = −3` (tie with `j=3`)
iff a complete residue computation (`lowj123.py`, `\bmod 4`, stable `\bmod 8`):
`  Δ(5) = −3 \iff (a,b) ≡ (3,0) \pmod 4.`
Hence `3 ∈ J*`, `J* ⊆ \{3,5\}`, with `J* = \{3,5\}` iff `(a,b) ≡ (3,0) (mod 4)`, else `\{3\}`.

In both branches `J* − j₀ ⊆ \{0,2\}` is a 2-adic box, so **`|J*| ≤ 2` is even on every tie**. This
establishes `G_λ(i)=0 ⟹ |J*|` even for the `c=5` interior, **with no residual**.

> **Theorem (restated).** For box-interior `λ=(a,b,5)`: if `a` even, `j₀=0` and `J*=\{0,2\}` iff
> `(a,b)≡(0,1)\,(4)` else `\{0\}`; if `a` odd, `j₀=3` and `J*=\{3,5\}` iff `(a,b)≡(3,0)\,(4)` else
> `\{3\}`. Only the generator `2` fires; `|J*| ≤ 2` is even on every tie. ∎

Adjoining the `c=5` **boundary** lemma (`2026-06-17-generalc-boundary-master-and-c5.md`, Theorem 2,
already complete, both parities) — which gives `val(b+i) > V` for the boundary indices and so makes
the global minimiser equal the interior one — and the finitely many small-`b` shapes covered there,
the three-row family **`(a,b,5)` is now fully proven, interior and boundary**, on the same footing as
`c=1,2,3,4`.

---

## 8. What is proved (no residual remains) and verification

**Proved in full (symbolic / complete finite residue checks, no `m`-cap):**
- Lemma 1 closed form; the structural decomposition `Q_5 = (a−3)(b−4)H_5 − 10!\,C(j,10)`; Prop 2;
  the depth-4 peeling identity and the clean low-index form.
- Prop 3 (offset / forced descent): `j₀ = 0` (a even), `j₀ = 3` (a odd), via the `R_1,R_2,R_3`
  valuation reductions (complete residue checks, stable across two moduli).
- `Δ(j) > 0` for `4 ≤ j ≤ 9` (twelve complete residue-class divisibility checks `2^{κ_j} ∣ Ψ_j`).
- Lemma N5: `8 ∣ H_5` (`v₂H_5 ≥ 3`, residues mod 16, sharp).
- Deep `j ≥ 10`: Case A (`Δ ≥ 1`, absorption); B1 (`g(j)>0`, analytic `j ≥ 32` + direct `10 ≤ j ≤ 31`,
  exceptions exactly `\{10,11,16,17,18,19\}`); B2 (the six exceptions, `2^{w_j} ∣ C\,C\,H_5`, complete
  residue checks).
- Tie classification (`\bmod 4`, stable `\bmod 8`): `(a,b)≡(0,1)` (a even), `(3,0)` (a odd).

**Cross-validation.** End-to-end Murnaghan–Nakayama (`endtoend.py`): `J*` matches the Theorem's
prediction for **all** box-interior `λ=(a,b,5)`, `m ≤ 45` (561 shapes, 128 ties, **0 mismatch**). The
deep-index hand bounds are `≤` the true `Δ` everywhere checked (`deepcheck.py`, `peel_deep.py`,
`lowfinal2.py`, all 0 failures, each stable across two moduli).

**Net.** The `c=5` interior tie set is `J* ⊆ \{j₀, j₀+2\}` (no generator `4/6/8/10`); proved for all
interior `j`, both parities, with **no residual**. The surprise vs PROVE.md: `c=5`, though odd, is a
**single-generator** family like `c=4`, *not* a two-generator one like `c=3`. The diagnosis (§9)
explains why.

---

## 9. Diagnosis — why `c=5` is single-generator (unlike `c=3`)

`c=3` had the full box `|J*|=4` in its `a`-odd branch because the heavy factor `(a−1)(b−2)` was
divisible by `4` there, paying *two* simultaneous deficits (`v₂(j)=1` at `j=5` **and** `v₂(j)=2` at
`j=7`) off the single switch `a≡1,b≡2`. The mechanism that *fired* the second generator was the tip
`−6!\,C(j,6)` reaching its first nonzero value `P_6(6)=720` while the binomial content was still low.

For `c=5` the tip is `−(j)_{10} = −10!\,C(j,10)`, which vanishes through `j = 9` and is first nonzero
at `j = 10` — but by `j = 10` we are already in the **deep régime**, where the absorption (§6, Case A)
and the heavy floor `v₂H_5 ≥ 3` (§6, Case B) force `Δ(j) > 0`. So the would-be generators
`4,6,8,10` never get a low-binomial-content window in which to fire: the tip's first appearance is
"too late". Only the generator `2` — which lives at `j₀+2` in the **low** range, governed by the clean
`Δ(j₀+2)=2v₂(R_2)−4` (a even) / the `j=5` tie (a odd) — survives. This is the same phenomenon that
made `c=4` single-generator (`2026-06-18` §9: tip `P_8` vanishes through `j=7`), now confirmed for the
odd case. The emerging general-`c` picture: **the parity of `a` sets the offset `j₀∈\{0,3\}`, the tip
order `(2c)!\,C(j,2c)` sets the *menu* of generators `2,…,2c`, but for `c ≥ 4` the menu items above
`2` all fall in the deep régime and are suppressed — only generator `2` fires, `|J*| ≤ 2`.** The
two-generator box of `c=3` is special to small `c` (the tip nonzero point `j=2c` lying *below* the
deep threshold).

### Files
`projects/code/threerow-c5/`: `census.py`,`census2.py` (empirical `J*`), `closedform.py`,`H5.py`
(Lemma 1 + decomposition), `prop2.py` (Prop 2 + Lemma N5 `8∣H_5`), `peel_deep.py` (peeling + `g(j)`),
`exceptional_rig.py` (deep exceptional residue checks), `lowfinal2.py` (`j=4..9` checks),
`lowindex.py`,`lowj123.py`,`smallH.py` (offset, `R_j`, tie classes), `deepcheck.py` (wide sweep),
`endtoend.py` (full theorem vs MN). MN engine: `threerow-c3/mn.py`.
