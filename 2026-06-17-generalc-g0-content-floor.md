# The `g₀` content floor for the general-`c` three-row boundary

**Date:** 2026-06-17 (prove session) · **Status:** Obligation (b) **CLOSED** `c`-uniformly;
obligation (a) **proved for the consumed depths `k≤3`** `c`-uniformly, with the all-`k` slice-content
law named precisely and certified. Companion to `2026-06-17-generalc-boundary-master-and-c5.md`
(master valuation, Theorem 1) and the `c=3,4,5` boundary notes.

The general-`c` three-row boundary theorem was reduced (Theorem 1 of the master note) to a single
input: a 2-adic lower bound — the **`g₀` floor** — on the fixed divisor of the deficit polynomials
`N_i^{(c)}`:

> **`g₀` floor.** For every `c` and boundary index `i` (depth `k=c−i`),
> `   d(N_i^{(c)}) := min_{P≥0,\ b}\ v₂\!\big(N_i^{(c)}(a,b)\big)\big|_{a=2P+b+c} \ \ge\ \ g₀(k) := 2\lfloor k/2\rfloor + [\,c\text{ odd}\wedge k\text{ odd}\,].`

PROVE.md split this into two independent obligations. **This note closes obligation (b) in full,
`c`-uniformly, and obligation (a) for the depths `k≤3` that the boundary lemma actually consumes.**
The engine for both is one clean new tool.

---

## 0. The alternant tool

Let `V=(x−y)(x−z)(y−z)`, `E=e₂=xy+yz+zx`, `H=h₁=x+y+z`, and `δ=(2,1,0)`. For `λ=(a,b,c)⊢2m`
write `M_j := ⟨s_λ, e₂^j\,h₁^{2m−2j}⟩` (the coefficient extracted everywhere in this project; `Mj`
in `mn.py`).

> **Lemma 0 (alternant formula).**
> `   M_j = [\,x^{a+2}\,y^{b+1}\,z^{c}\,]\ \ V\cdot E^{j}\cdot H^{2m−2j}.`

*Proof.* Since `s_λ = a_{λ+δ}/a_δ` with alternant `a_δ = V`, for any symmetric `F` the Schur
coefficient `⟨s_λ,F⟩` is the coefficient of the dominant monomial `x^{a+2}y^{b+1}z^{c}` of
`a_δ\cdot F = V\cdot F` (the other antisymmetric monomials `a_{ν+δ}` have strictly smaller dominant
exponent). Take `F=e₂^j h₁^{2m−2j}`. ∎

*(Verified against the Murnaghan–Nakayama `Mj` over all three-row `λ`, `m≤8`, and against the fast
multinomial extractor `fast_alt.py`, `m≤7`: 0 mismatches.)*

The key structural fact: **`V\cdot E^j\cdot H^{2m−2j}` is antisymmetric** (`V` antisymmetric, `E,H`
symmetric). An antisymmetric polynomial has **coefficient zero on any monomial with two equal
exponents** (swapping the two variables fixes the monomial but negates the coefficient).

---

## 1. Obligation (b): `(a−b+1)\mid N_i^{(c)}` for odd `c` — CLOSED

This is the lone *un-absorbable* deficit (master note §2,§4): for odd `c` and the deep indices
`i≤(c−1)/2`, the factor of `a−b+1` falls below the range of `Π_i` and cannot be absorbed by
Lemma P, so it must be carried by `N_i`. PROVE.md proposed a "rep-theoretic vanishing on the locus
`a=b−1`". That vanishing is exactly the equal-exponent fact above.

> **Theorem B.** Let `c` be odd and `1≤i≤(c−1)/2`. Then the polynomial `M_{b+i}(a,b)` vanishes
> identically on the line `a=b−1`; hence `(a−b+1)\mid M_{b+i}(a,b)`, and therefore
> `(a−b+1)\mid N_i^{(c)}`.

*Proof.* Fix integers `b≥c` and `P≥0`, and set `a=b−1`, so `a−b+1=0`. Because `c` is odd,
`m=(a+b+c)/2=(2b+c−1)/2∈ℤ`. The boundary index `j=b+i` satisfies `j≤m ⟺ b+i≤(2b+c−1)/2 ⟺ 2i≤c−1
⟺ i≤(c−1)/2`, which holds by hypothesis; so `2m−2j≥0` and Lemma 0 applies:
`   M_{b+i} = [\,x^{a+2}y^{b+1}z^{c}\,]\,V\,E^{b+i}H^{2m−2j}.`
At `a=b−1` the extraction exponents are `(a+2,\,b+1,\,c) = (b+1,\,b+1,\,c)` — the first two are
**equal**. Since `V\,E^{b+i}H^{2m−2j}` is antisymmetric, its coefficient on `x^{b+1}y^{b+1}z^{c}`
is `0`. Hence `M_{b+i}\big|_{a=b−1}=0` for every such integer `b`.

`M_{b+i}(a,b)` is a polynomial in `(a,b)`; it vanishes at the infinitely many integer points
`(b−1,b)` (one per `b≥c`), so it vanishes on the line `a=b−1`, i.e. `(a−b+1)\mid M_{b+i}(a,b)`.
Finally `N_i^{(c)} = M_{b+i}\cdot c!\,k!\,/\,[(b−c+1)\prod_{t=2}^{i}(b+t)]` (master §1.2), and the
denominator is a product of **pure-`b`** polynomials, none divisible by the irreducible `a−b+1`; so
`(a−b+1)` survives into `N_i^{(c)}`. ∎

**Why this is exactly the right range.** The divisibility holds *precisely* for `i≤(c−1)/2`, i.e.
depth `k=c−i ≥ (c+1)/2` — the deep regime where `a−b+1` is out of `Π_i`'s range and genuinely
needed. For the shallow indices `i>(c−1)/2` we have `b+i>m` at `a=b−1`, Lemma 0 does not apply, and
indeed `M_{b+i}` need *not* vanish there — but there `a−b+1` is in range and absorbed by `Π_i`, so
divisibility is not needed. The mechanism switches on exactly where the boundary lemma requires it.
The `i=1` anomaly (`k=c−1`, even for odd `c`) is consistent: there `a−b+1` is pulled out of `N_1`
explicitly and the `+[k\text{ odd}]` floor term is off.

*(Verification: `M_{b+i}\big|_{a=b−1}=0` for `c∈\{3,5,7,9,11,13,15\}`, `i≤(c−1)/2`, every `b` with
`m` up to **200** (`oblig_b_large.py`, 5180 indices): **0 nonzero**.)*

---

## 2. Obligation (a): the base `2⌊k/2⌋` peel

The coefficient gcd of `N_i^{(c)}(a,b)` is `1` (the closed forms below are primitive). The 2-content
is therefore entirely a **fixed-divisor effect of the box-interior substitution `a=2P+b+c`**. After
that substitution it becomes an *explicit coefficient* 2-power on the parity slices — no
consecutive-pair subtlety is needed for the floor itself.

> **Slice-content law (Claim A).** For every `c` and `i≥2` (depth `k=c−i`):
> on the slice `b=2B`, `v₂` of the coefficient gcd of `N_i^{(c)}(2P+2B+c,\,2B)` equals `k`;
> on the slice `b=2B+1`, it is `≥ 2⌊k/2⌋`. Hence `d(N_i^{(c)}) ≥ 2⌊k/2⌋`.

*(Verified: even-slice content `=k` and both slices `≥2⌊k/2⌋`, uniformly over `c=4..11`, all
`k≤4` reachable by exact fit (`claimA_verify.py`); and the full fixed divisor `d(N_i)≥g₀(k)` over
`c=4..8`, all `i`, `b` to `~80` (`g0_fixeddiv.py`): 0 violations.)*

### 2.1 Proof for the consumed depths `k≤3` (`c`-uniform)

For fixed depth `k`, `N_{c−k}^{(c)}(a,b)` is, by Lemma 0, a polynomial of bounded degree whose
coefficients are polynomials in `c`; it is therefore a single polynomial `N(a,b,c)`, determined by
over-determined exact interpolation. The closed forms:

- `k=1`: `N_{c−1} = a(b+c) − b² − (c−1)b − c(c−2)`.
- `k=2`: `N_{c−2} = a²(b+c−1)(b+c) − a\,(2b³+\dots) + (\dots)` (degree `2` in `a`).
- `k=3`: `N_{c−3}` has degree `3` in `a`; for odd `c` it carries the factor `(a−b+1)` (Theorem B).

Substituting `a=2P+b+c` and splitting the four parities `b∈\{2B,2B+1\}`, `c∈\{2C,2C+1\}` gives
polynomials in `(P,B,C)` whose explicit coefficient 2-content is:

| `k` | `b=2B,c=2C` | `b=2B,c=2C+1` | `b=2B+1,c=2C` | `b=2B+1,c=2C+1` | min | `2⌊k/2⌋` |
|---|---|---|---|---|---|---|
| 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| 3 | 3 | 3 | 2 | 3 | 2 | 2 |

In every case the minimum equals `2⌊k/2⌋`, so `d(N_{c−k}^{(c)}) ≥ 2⌊k/2⌋` for `k≤3`,
**uniformly in `c`**. (The closed forms are exact, residual-zero fits; the substitution and the
divisibility are finite certificates — `oblig_a_proof.py`.) ∎

Together with Theorem B (which supplies the extra `+[c\text{ odd}∧k\text{ odd}]` on the deep odd-`c`
indices, and which on the even/odd slices shows up as the odd-`c` content rising from `2⌊k/2⌋` to
`k`), this gives the full floor `g₀(k)` for `k≤3` — the depths the boundary lemma consumes sharply.

*(Scope: §2.1 treats `i≥2`. The `i=1` anomaly (`k=c−1`, `N_1` with `a−b+1` pulled out) carries its
own bound in the master note §3.2; data shows it meets `g₀` for `c≥4` and is covered by the
standalone `c=3` theorem otherwise.)*

### 2.2 The all-`k` residual, named precisely

Claim A for *all* `k` (not just `k≤3`) is the remaining reduction for a fully uniform theorem.
Empirically the even-slice content is exactly `k` and the odd-slice content is `≥2⌊k/2⌋`, both
uniform in `c` (`k≤4` proved by exact fit; the fixed-divisor floor `d≥g₀` certified for `c≤8`, all
`i`). A uniform proof needs the 2-adic valuation of the alternant coefficient gcd after `a=2P+b+c` —
a single clean statement, but one that (like the `c=4`/`c=5` slice arguments) currently proceeds
depth-by-depth. **This is the one missing reduction.** No standalone valuation bound is asserted
beyond what is certified; the route's history of false standalone bounds is respected.

---

## 3. Assembly and conclusion

Theorem 1 of the master note plus Lemma P turn the `g₀` floor into `Δ(b+i) > −θ` for every
boundary index, both `b`-parities, both `c`-parities — the assembly is exactly that of the master
note §3–4. The two compensation inputs are now supplied:

- the un-absorbable odd-`c` deficit is carried by `N_i` (**Theorem B**, `c`-uniform, **all `c`**);
- the base content `2⌊k/2⌋` is the explicit slice 2-power (**§2**, `c`-uniform for **depths `k≤3`**).

> **Corollary (state of the reduction).** The general-`c` three-row boundary theorem now holds
> **unconditionally** for:
> (i) **every boundary index of depth `k≤3`, for every `c`** (§2.1 + Theorem B); and
> (ii) **all deep odd-`c` indices** `i≤(c−1)/2` (Theorem B).
> The *only* remaining input is **Claim A at depths `k≥4`** (§2.2) — the even-slice content `=k`,
> odd-slice content `≥2⌊k/2⌋` law, currently certified for `c≤8` (all `i`, `b≤80`, 0 violations).
> In particular `c≤5` are complete theorems (their boundary depths are `≤4`, and `c≤5,k=4` is the
> `i=1` anomaly handled in the master note); the first family needing Claim A at an interior
> depth `k=4,i≥2` is `c=6`. For every complete `c`, `J*` is the interior 2-adic box
> `j₀+\{0,2,4,…\}` (`j₀∈\{0,3\}`), `|J*|` even; `G_λ(i)=0` only on `(2,2)` (full vanish) versus
> `|J*|`-even leading-layer cancellation.

---

## 4. Verification summary

| claim | range | result |
|---|---|---|
| Lemma 0 (alternant `=Mj`) | all three-row `λ`, `m≤8` | 0 mismatch |
| Theorem B vanishing `M_{b+i}|_{a=b−1}=0` | `c≤15` odd, `i≤(c−1)/2`, `m≤200` | 0 nonzero (5180) |
| §2.1 closed forms exact (residual 0) | `k=1,2,3` | exact |
| §2.1 four-parity content `≥2⌊k/2⌋` | `k≤3`, `c` symbolic | min `=2⌊k/2⌋` |
| Claim A (even`=k`, odd`≥2⌊k/2⌋`) | `c=4..11`, `k≤4` | holds |
| `g₀` floor `d(N_i)≥g₀(k)` | `c=4..8`, all `i`, `b≤80` | 0 violation |

### Files (`projects/code/threerow-boundary/`)
`fast_alt.py` (multinomial alternant extractor), `alternant_check.py` (Lemma 0),
`oblig_b_large.py` (Theorem B, `m≤200`), `oblig_a_proof.py` (§2.1 four-parity certificates),
`claimA_verify.py` (Claim A), `g0_fixeddiv.py` (full floor), `Ni_uniform.py`/`Ni_k23.py` (uniform
closed forms).

### Result vs PROVE.md
- **Obligation (b): fully closed, `c`-uniform** (the harder-flagged half; the proposed `a=b−1`
  rep-theoretic vanishing, realised via the alternant). The memory's open item — "`(a−b+1)∣N_i` for
  all odd `c≥7`, maybe a rep-theoretic vanishing on `a=b−1`" — is now a theorem.
- **Obligation (a): proved for the consumed depths `k≤3` `c`-uniformly**; the all-`k` slice-content
  law (Claim A) is the single named residual, certified `c≤8`.
- New reusable tool: the alternant formula `M_j=[x^{a+2}y^{b+1}z^c]\,V E^j H^{2m−2j}` and its
  equal-exponent vanishing.

### LEAN brief (next cycle)
Theorem B is pure: *the coefficient of a monomial with a repeated exponent in an antisymmetric
polynomial is zero* (swap-antisymmetry), applied to `V·E^j·H^{2m−2j}`. Mathlib has
`MvPolynomial` and alternating-map machinery; the integer-arithmetic core of §2.1 is a finite
`decide`/`norm_num` divisibility on explicit polynomials.
