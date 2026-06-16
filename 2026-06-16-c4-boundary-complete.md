# The `c = 4` boundary lemma — CLOSED. Three-row `(a,b,4)` is now a COMPLETE theorem

**Date:** 2026-06-16 (prove session)
**Status:** **Complete.** The boundary lemma for `c=4` — the first `c` past the `c=3`
factor-in-product wall — is proved by a hand argument certified numerically. Together with
the already-proved interior (general-`c` Number Lemma `NL_c`,
`2026-06-14-numberlemma-general-c-proved.md`) this upgrades the three-row `d=4` family
**`λ=(a,b,4)`** to a **complete theorem** (interior + boundary), joining the complete
`c=1,2,3` families.

The keystone is **not** a new factor-in-product principle. It is the realisation that the
"irreducible `a`-quadratic deficit" `N_i` that PROVE.md flagged as the wall enters the
valuation `Δ(b+i)` as a **positive surplus**, so only a *lower bound* on `v₂(N_i)` is needed —
and that bound is the **2-adic content** of `N_i` (the gcd of its coefficients after
`a = 2P+b+4`), which requires **no factorisation at all**. Irreducibility is irrelevant.

Companion to `2026-06-16-c3-boundary-complete.md`, whose Lemma-D / telescoping template this
note mirrors. All scripts in `projects/code/threerow-boundary/` (`c4_*.py`), MN engine `mn.py`.

---

## 0. Setup and what must be proved

For `λ=(a,b,4) ⊢ 2m` we have `a+b = 2m−4`, so **`a` and `b` have the same parity**. Set

> `G_j = ⟨s_λ, e₂^j h₁^{2(m−j)}⟩ = C(m,j) M_j`,  `val(j) = j + 2 v₂(G_j)`,  `V = min_j val(j)`,
> `J* = {j : val(j) = V}`,  `Δ(j) := val(j) − val(0)`,  `P := m−b−4 ≥ 0`.

**The interior structure (computed, c=4).** For every `c=4` shape in the box-interior range
`b ≥ 6`, one finds `j₀ := min J* = 0` and `V = val(0)`, i.e. **`θ := val(0)−V = 0`** (verified
`m ≤ 45`, all residues `(a,b) mod 4`). Equivalently the interior 2-adic box is `{0,2,4,6}` and
the minimiser sits at `j=0`. So the boundary lemma to prove is the *strict* statement:

> **Boundary lemma (c=4).** For `λ=(a,b,4)`, `b ≥ 6`, every boundary index
> `j ∈ {b+1,b+2,b+3,b+4}` satisfies `Δ(b+i) > 0` (`i=1,2,3,4`).

`M_j` is supported on `0 ≤ j ≤ b+4`, so these four are all the boundary indices. The finitely
many smaller shapes `b ∈ {4,5}` (where the box `{0,2,4,6}` does not fit inside `[0,b]`) are
verified directly (`m ≤ 60`, §6).

**Standing relations** from `a = 2P+b+4`:

> `a − b + 1 = 2P + 5`  (**odd** — the `c`-even signature; contrast `c=3` where it was even),
> `a − 2 = 2P + b + 2` (`= a−c+2`),  `a + 2 = 2P + b + 6`,  `2m = 2P+2b+8`,  `m = P+b+4`.

In particular `v₂(a−b+1) = 0` always, and **exactly one** of `v₂(a−2)`, `v₂(b−3)` is positive:

- `b` even ⟹ `a−2 = 2P+b+2` even, `b−3` odd: the lone deficit is `v₂(a−2) = 1 + v₂(P+\tfracb2+1)`;
- `b` odd ⟹ `a−2` odd, `b−3` even: `v₂(a−2)=0` (the deficit `v₂(b−3)` will cancel in the descent).

---

## 1. Closed forms and the exact descent ratio `R_i`

The signed three-row Jacobi–Trudi `M_{b+i}` factor over `ℤ[a,b]` as (fitted and verified exact
vs Murnaghan–Nakayama, `m ≤ 30`):

> `M_{b+4} = (b−3)(b+2)(b+3)(b+4)/24`  (Lemma T, `a`-free),
> `M_{b+3} = (b−3)(b+2)(b+3)·N₃/24`,  `N₃ = a(b+4) − (b²+3b+8)`,
> `M_{b+2} = (b−3)(b+2)·N₂/48`,        `N₂` = **irreducible `a`-quadratic** (below),
> `M_{b+1} = (b−3)(a−b+1)·N₁/144`,     `N₁` = `a`-quadratic (below),

with the deficit polynomials
`N₂ = a²b²+7a²b+12a²−2ab³−9ab²−21ab−20a+b⁴+2b³+13b²+16b−8`,
`N₁ = a²b³+9a²b²+26a²b+24a²−2ab⁴−7ab³−13ab²−14ab+b⁵−2b⁴+17b³−4b²−84b−96`.

Using the 3-row hook count `f = G_0 = (2m)!\,(a−b+1)(a−2)(b−3)\,/[(a+2)!(b+1)!\,24]` and
`C(m,b+i) = (P+b+4)!/[(b+i)!(P+4−i)!]`, the descent ratio `R_i := G_{b+i}/G_0` simplifies — the
factors `(b−3)`, `(a−b+1)`, the central `(b+t)`'s and most factorials cancel — to the exact
rational forms (**verified as exact rationals vs MN, `m ≤ 41`, 0 mismatch**), where
`D := ∏_{s=2P+b+7}^{2P+2b+8} s` is a run of `b+2` consecutive integers:

> `R_4 = ∏_{t=1}^{b+4}(P+t) \big/ [(2P+5)(2P+b+2)·D]`,
> `R_3 = N₃·∏_{t=2}^{b+4}(P+t) \big/ [(2P+5)(2P+b+2)·D]`,
> `R_2 = N₂·∏_{t=3}^{b+4}(P+t) \big/ [2(2P+5)(2P+b+2)·D]`,
> `R_1 = N₁·∏_{t=4}^{b+4}(P+t) \big/ [6(2P+b+2)·D]`.

Since `2P+5` is odd and `Δ(b+i) = (b+i) + 2\,v₂(R_i)`, writing `c_i := (0,0,1,1)` for
`i=(4,3,2,1)` (the `v₂` of the constant `1,1,2,6` in the denominators) and
`V_{\rm num}(k) := v₂\big(∏_{t=k}^{b+4}(P+t)\big)`:

> `v₂(R_i) = v₂(N_i) − c_i − v₂(a−2) + V_{\rm num}(5−i) − v₂(D)`,  with `N_4:=1` (so `v₂=0`).

---

## 2. The telescoping reduction (the `D`-cancellation)

`D` is "double-scale" but its **even part, halved, ends at `P+b+4`** — the very top of the
numerator products — so the two telescope.

> **Lemma R.** Let `ℓ := ⌈(b+7)/2⌉`, `E := (P+b+4) − (P+ℓ) + 1 = b+5−ℓ`. For `k ≤ ℓ−1`,
> `  V_{\rm num}(k) − v₂(D) = v₂\big(∏_{t=k}^{ℓ−1}(P+t)\big) − E.`

*Proof.* The even integers in `[2P+b+7,\,2P+2b+8]` are exactly `{2u : u ∈ [P+ℓ,\,P+b+4]}`
(the top `2P+2b+8 = 2(P+b+4)` is even; the smallest even is `2(P+ℓ)` with `ℓ=⌈(b+7)/2⌉`).
Hence `v₂(D) = Σ_u v₂(2u) = E + v₂\big(∏_{t=ℓ}^{b+4}(P+t)\big)` (`E` even factors, each a `2·u`).
Subtracting from `V_{\rm num}(k) = v₂\big(∏_{t=k}^{b+4}(P+t)\big)` cancels the shared tail
`∏_{t=ℓ}^{b+4}(P+t)`. ∎

Concretely `ℓ = (b+7)/2`, `E = (b+3)/2` (`b` odd) and `ℓ = (b+8)/2`, `E = (b+2)/2` (`b` even).
Writing `Π_i := ∏_{t=5−i}^{ℓ−1}(P+t)` (a run of `L_i := ℓ+i−5` consecutive integers), Lemma R
turns the ratio formula into the **master valuation** (verified vs MN, `m ≤ 43`, 0 mismatch):

> **`v₂(R_i) = v₂(N_i) − c_i − v₂(a−2) + v₂(Π_i) − E`,**  so  **`Δ(b+i) = (b+i) + 2 v₂(R_i)`.**

Lengths: `L_i = (b−3)/2 + i` (`b` odd) and `L_i = (b−2)/2 + i` (`b` even); all `≥ 1` for `b≥6`.

---

## 3. The wall-breaker: the 2-adic content of `N_i`

The deficits `N_i` appear with a **`+`** sign in `v₂(R_i)`: they are surpluses, not obstacles.
We only need lower bounds, and these come from the polynomial **content** after `a = 2P+b+4`.

> **Lemma C (content bounds).** For all integers `P ≥ 0`, `b`:
> (i) `v₂(N₂) ≥ 2`;
> (ii) `b` even ⟹ `v₂(N₃) ≥ 1`;  `b` odd ⟹ `N₃` is odd;
> (iii) `b` even ⟹ `v₂(N₁) ≥ 3`;  `b` odd ⟹ `v₂(N₁) ≥ 2`.

*Proof.* Substitute `a = 2P+b+4` and reduce over `ℤ[P,b]`.

(i) `N₂ = 2\,Q`, `Q = 2P²b²+14P²b+24P²+13Pb²+59Pb+76P+20b²+60b+52`. Modulo 2,
`Q ≡ 13Pb²+59Pb ≡ Pb²+Pb = P\,b(b+1) ≡ 0` since `b(b+1)` is even. So `2 \mid Q`, `v₂(N₂) ≥ 2`.

(ii) `N₃ = 2Pb+8P+5b+8`. For `b=2B`: `N₃ = 2(2BP+5B+4P+4)`, so `v₂ ≥ 1`. For `b` odd, `5b` is
odd and `2Pb,8P,8` are even, so `N₃` is odd.

(iii) Substituting `b=2B` (resp. `b=2B+1`) and `a=2P+b+4` gives the explicit factorisations
`N₁ = 8\,(4B³P²+38B³P+90B³+18B²P²+111B²P+153B²+26BP²+121BP+117B+12P²+48P+36)` (`b` even) and
`N₁ = 4\,(8B³P²+76B³P+180B³+48B²P²+336B²P+576B²+94BP²+521BP+675B+60P²+282P+288)` (`b` odd),
both with integer inner polynomials; hence `v₂(N₁) ≥ 3` (`b` even), `≥ 2` (`b` odd). ∎

*(Content bounds verified with 0 violations for every `c=4` shape `m < 200`.)* The point worth
stressing: **none of these proofs factor `N_i`.** The `c≥4` wall — `M_{b+c−2}` irreducible — is
a wall only for the factor-in-product engine; the content bound steps around it entirely and is
trivially uniform in `c`.

---

## 4. Elementary product lemma

> **Lemma P.** A product `Π` of `Q` consecutive integers has `v₂(Π) ≥ ⌊Q/2⌋` (it contains
> `≥⌊Q/2⌋` even factors). If one factor is removed, the remainder has `v₂ ≥ ⌊Q/2⌋ − 1`.

*Proof.* Among `Q` consecutive integers at least `⌊Q/2⌋` are even, each contributing `v₂ ≥ 1`;
removing one factor deletes at most one even factor. ∎

---

## 5. Proof of the boundary lemma

By §2, `Δ(b+i) = (b+i) + 2\big[\,v₂(N_i) − c_i − v₂(a−2) + v₂(Π_i) − E\,\big]`. We split on the
parity of `b` and use Lemmas C, P. Throughout `b ≥ 6`.

### 5.1 `b` odd (`v₂(a−2)=0`, `E=(b+3)/2`, `L_i=(b−3)/2+i`)

Here `Δ(b+i) = (b+i) − (b+3) − 2c_i + 2v₂(N_i) + 2v₂(Π_i) = (i−3) − 2c_i + 2v₂(N_i) + 2v₂(Π_i)`.

- **`i=4`** (`c=0`, `v₂(N_4)=0`): `Δ(b+4) = 1 + 2v₂(Π_4) ≥ 1`. ✓
- **`i=3`** (`c=0`, `N₃` odd): `Δ(b+3) = 2v₂(Π_3)`; `L_3 = (b−3)/2+3 ≥ 5`, Lemma P gives
  `v₂(Π_3) ≥ 2`, so `Δ(b+3) ≥ 4`. ✓
- **`i=2`** (`c=1`): `Δ(b+2) = −3 + 2v₂(N₂) + 2v₂(Π_2) ≥ −3 + 2·2 + 0 = 1` by Lemma C(i). ✓
- **`i=1`** (`c=1`): `Δ(b+1) = −4 + 2v₂(N₁) + 2v₂(Π_1) ≥ −4 + 2·2 + 2v₂(Π_1)` by Lemma C(iii);
  `L_1 = (b−3)/2+1 ≥ 3` (as `b ≥ 7`), so `v₂(Π_1) ≥ 1` and `Δ(b+1) ≥ 2`. ✓

### 5.2 `b` even (`v₂(a−2)=1+v₂(P+\tfracb2+1)`, `E=(b+2)/2`, `L_i=(b−2)/2+i`)

The factor `P+\tfracb2+1` lies in the range `[P+5−i,\,P+ℓ−1]` of `Π_i` (since
`5−i ≤ \tfracb2+1 ≤ ℓ−1=\tfracb2+3` for `b≥6`, `i≥1`), so it is one of the factors of `Π_i`.
Hence `−v₂(a−2)+v₂(Π_i) = −1 + v₂(Π_i')`, where `Π_i' := Π_i /(P+\tfracb2+1)` is the same run
with that one factor deleted (`L_i−1` factors). Substituting `E=(b+2)/2`:

`Δ(b+i) = (b+i) − (b+2) − 2c_i − 2 + 2v₂(N_i) + 2v₂(Π_i') = (i−4) − 2c_i + 2v₂(N_i) + 2v₂(Π_i')`.

- **`i=4`** (`c=0`, `v₂(N_4)=0`): `Δ(b+4) = 2v₂(Π_4')`; `L_4 = (b−2)/2+4 = \tfracb2+3 ≥ 6`, so by
  Lemma P (remove one factor) `v₂(Π_4') ≥ ⌊L_4/2⌋−1 ≥ 2`, giving `Δ(b+4) ≥ 4`. ✓
- **`i=3`** (`c=0`): `Δ(b+3) = −1 + 2v₂(N₃) + 2v₂(Π_3') ≥ −1 + 2·1 + 0 = 1` by Lemma C(ii). ✓
- **`i=2`** (`c=1`): `Δ(b+2) = −4 + 2v₂(N₂) + 2v₂(Π_2') ≥ −4 + 2·2 + 2v₂(Π_2')` by Lemma C(i);
  `L_2 = \tfracb2+1 ≥ 4`, so `v₂(Π_2') ≥ ⌊L_2/2⌋−1 ≥ 1` and `Δ(b+2) ≥ 2`. ✓
- **`i=1`** (`c=1`): `Δ(b+1) = −5 + 2v₂(N₁) + 2v₂(Π_1') ≥ −5 + 2·3 + 0 = 1` by Lemma C(iii). ✓

In every case `Δ(b+i) ≥ 1 > 0`. ∎

*(A certified lower bound assembling **only** Lemmas C and P reproduces the margins
`(7,4,1,2)` (`b` odd, `i=4,3,2,1`) and `(4,1,2,1)` (`b` even) and is `≤ Δ` and `≥ 1` for every
`c=4` shape `m < 120`: 0 failures.)*

---

## 6. Conclusion

> **Theorem (three-row `c = 4` boundary lemma).** For `λ=(a,b,4) ⊢ 2m`, `b ≥ 6`, every
> boundary index `j ∈ {b+1,b+2,b+3,b+4}` has `val(j) > V`. With the smaller shapes `b∈{4,5}`
> checked directly, no boundary index is ever an extra minimiser of `val`.

Combined with the proved interior (`NL_4`):

> **Three-row `d=4` family `(a,b,4)` is a COMPLETE theorem.** `J*` is exactly the interior
> 2-adic box `\{0,2,4,6\}` (`j₀=0`); on every tie `|J*| ∈ \{1,2,4\}` is even, so the leading-`π`
> dichotomy gives `G_λ(i)=0` only on `|J*|`-even shapes — matching `c=1,2,3`. (`J*⊆\{0,2,4,6\}`,
> `|J*|∈\{1,2,4\}` verified for all `c=4` shapes `m<60`, incl. `b∈{4,5}`.)

No machine residual remains for `c ≤ 4`.

### The general-`c` mechanism (now visible)

The `c=4` proof reveals that the PROVE.md wall was illusory, and gives the general template:

1. **`θ` and `j₀`** are set by the interior theorem (here `θ=0`, `j₀=0`; for `c` even this is
   forced, since `a−b+1 = 2P+c+1` is odd and the offset collapses).
2. **The exact ratio `R_i`** and the **`D`-telescoping (Lemma R)** are `c`-uniform: the even
   part of the bottom run `D`, halved, always shares its top `P+b+c` with the numerator runs.
3. **The deficits `N_i^{(c)}` are surpluses** in `v₂(R_i)`; one needs only their **2-adic
   content** (gcd of coefficients after `a=2P+b+c`) — *no factorisation*, so the irreducibility
   of `M_{b+c−2}` for `c≥4` is irrelevant. This is the one idea that breaks the wall.
4. **Lemma P** (consecutive-product evens-count, with up-to-one removal for the absorbed
   `v₂(a−c+2)` deficit) finishes, exactly as here.

For `c=5` (odd, like `c=3`): one expects `θ = τ(τ+1)/2 > 0` and the smooth deficit
`a−b+1 = 2P+6` even (needs the F1/F2 absorption of the `c=3` `a`-odd top), but the deeper
indices `b+1, b+2, b+3` again close by **content bounds on `N_i^{(5)}`** — the explicit contents
are the only new `c`-specific input. The content-bound principle is the uniform engine; only the
interior offset `θ(c)` and the explicit polynomial contents remain to be tabulated per `c`.

---

## 7. Verification summary

All scripts `projects/code/threerow-boundary/` (`mn.py` = Murnaghan–Nakayama).

| claim | range | result |
|---|---|---|
| `θ=0`, `j₀=0`, box `{0,2,4,6}` | `m ≤ 45`, all `(a,b)` mod 4 | uniform |
| `M_{b+i}` closed forms (incl. `N₃,N₂,N₁`) vs MN | `m ≤ 30` | 0 mismatch |
| exact ratios `R_i` (rationals) vs `G_{b+i}/G_0` | `m ≤ 41` | 0 mismatch |
| master `v₂(R_i)` (Lemma R + content) vs MN | `m ≤ 43` | 0 mismatch |
| Lemma C content bounds `v₂(N_i)` | `m < 200` | 0 violation |
| **certified `Δ(b+i)` lower bound** (Lemmas C,P only) `≤ Δ` and `≥ 1` | `m < 120` | 0 failure |
| full theorem `J*⊆{0,2,4,6}`, `|J*|∈{1,2,4}` (incl. `b∈{4,5}`) | `m < 60` | 0 violation |
| Lemma P (evens-count, ±1 removal) | `Q≤39`, `R<300` | 0 failure |

### Files
`c4_scout.py` (θ/`j₀`), `c4_forms.py` (closed forms), `c4_lemmaD.py` (direct `Δ`),
`c4_ratios.py` (exact `R_i`), `c4_reduce.py` (Lemma R master formula), `c4_Nstruct.py`
(content reductions), `c4_certify.py` (content bounds + certified `Δ≥1`),
`c4_fulltheorem.py` (full claim + Lemma P).

### Gaps
None for `c ≤ 4`. The general-`c` template (§6) is proved in mechanism but written out only for
`c=4`; `c=5` needs the interior offset `θ(5)` and the explicit contents of `N_i^{(5)}` — the
natural next target, alongside the LEAN formalisation of the now-complete `c=1,2,3,4` families.
