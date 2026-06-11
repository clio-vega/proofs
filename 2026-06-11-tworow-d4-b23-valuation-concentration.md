# Two-row d=4 law, classes b ≡ 2,3 (mod 4): the 2-adic valuation concentrates in a single root

**Date:** 2026-06-11 (prove session)
**Status:** **NOT a uniform closure.** New, rigorous, and clean: an exact closed form for the
2-adic valuation `v₂(I_b(m))` valid at *every* integer `m`, which proves that the entire
(unbounded) 2-adic behaviour of `I_b` is carried by **one** term — the 2-adic distance
`v₂(m − ρ_b)` to the unique 2-adic root. This is the load-bearing lemma flagged in the
session brief; it reframes "v₂ unbounded, no truncation works" into "v₂ is one controllable
distance, and the whole law is the single statement `ρ_b ∉ ℤ_{≥b}`." Verified for all
`b ≡ 2,3 (mod 4)`, `6 ≤ b ≤ 80`, across thousands of integers `m`.

---

## 0. The object and what was already known

For `λ = (2m−b, b)`, `s = 1+i`, `A = 1+su+u²`, set
`G_b(m) = [u^b]((1−u)A^m)` and `I_b(m) = Im G_b(m)`. The two-row d=4 fibre law
`G_{(2m−b,b)}(i) = 0 ⟺ (2,2)` is equivalent to `I_b(m) ≠ 0` for every integer `m ≥ b`.
Classes `b ≡ 0,1 (mod 4)` are proved (06-07/06-08). This note concerns `b ≡ 2,3`.

Established in prior cycles (cited, not re-derived here):

- **Reduction.** `I_b(m) = (m)_{R+1} · Q_b(m)`, where `(m)_{R+1} = m(m−1)⋯(m−R)`,
  `R = ⌊(b−1)/2⌋`, and `Q_b ∈ ℚ[m]` has degree `d = ⌊b/2⌋`. The law for `b` is equivalent to
  `(★)`: `Q_b` has no integer root `m ≥ b`. (Note `I_b` has rational, integer-valued
  coefficients; so does `Q_b`. Its leading coefficient `ℓ_b ∈ ℚ` need not be a 2-adic integer.)
- **Lemma 1 (mod-2 collapse).** `I_b(m) ≡ m·C(m−1, R) (mod 2)`.
- **Lemma 2 (monic mod-2 factorisation).** The *monic* polynomial `Q̂_b := Q_b/ℓ_b` lies in
  `ℤ₂[X]` and reduces mod 2 to
  `Q̂_b ≡ X·(X²+X+1)^{(d−1)/2} (mod 2)`, with `X²+X+1` irreducible over `𝔽₂`
  (and `(d−1)/2 = ⌊(b−1)/4⌋ ∈ ℤ` because `d` is odd for `b ≡ 2,3`).
  *Re-verified this session for all `b ≡ 2,3`, `6 ≤ b ≤ 80`.* (The earlier worry that this
  "breaks" at `b = 26,34,…` was an artifact of reducing a non-monic primitive integer model;
  the monic model is clean throughout.)
- **Theorem 2 (unique 2-adic root).** By Lemma 2 + Hensel, `Q̂_b` (hence `Q_b`) has exactly one
  root `ρ_b ∈ ℤ₂`; it is simple and `v₂(ρ_b) = 1`. Consequently `Q_b` has at most one rational
  root, namely `ρ_b`, and `(★)` ⟺ `ρ_b ∉ ℤ_{≥b}`.

The sharp open gap is to prove `ρ_b ∉ ℤ_{≥b}` **uniformly in b**.

## 1. Lemma A — the valuation concentrates in a single distance

> **Lemma A.** Let `b ≡ 2,3 (mod 4)`, `R = ⌊(b−1)/2⌋`, and let `ρ_b ∈ ℤ₂` be the unique 2-adic
> root of `Q_b`. Then for **every** integer `m`,
> ```
>     v₂( I_b(m) )  =  v₂(ℓ_b)  +  Σ_{r=0}^{R} v₂(m − r)  +  v₂(m − ρ_b),
> ```
> where `ℓ_b ∈ ℚ` is the leading coefficient of `Q_b`. In particular the second term is the
> elementary Legendre/Kummer sum `v₂((m)_{R+1})`, bounded on each residue window, and the
> **entire unbounded part** of `v₂(I_b(m))` is the single 2-adic distance `v₂(m − ρ_b)`.

**Proof.** Factor `Q_b` over the algebraic closure `ℂ₂ = \overline{ℚ_2}`:
`Q_b(X) = ℓ_b ∏_{θ}(X − θ)`, the product over the `d` roots `θ`. Group the roots as `ρ_b`
together with the remaining `d−1` roots, which I call *non-special*. For any integer `m`,
`v₂(I_b(m)) = v₂((m)_{R+1}) + v₂(Q_b(m)) = Σ_{r=0}^R v₂(m−r) + v₂(ℓ_b) + Σ_θ v₂(m − θ)`.
It remains to show each non-special root contributes `v₂(m − θ) = 0` for every integer `m`.

By Lemma 2 the monic `Q̂_b` lies in `ℤ₂[X]` and reduces mod 2 to `X·(X²+X+1)^{(d−1)/2}`. Hence
every root `θ` of `Q̂_b` is integral over `ℤ₂` (a root of a monic `ℤ₂`-polynomial), so
`v₂(θ) ≥ 0`, and its reduction `\barθ` in the residue field `𝔽̄₂` is a root of
`X(X²+X+1)^{(d−1)/2}`. For a **non-special** root, `\barθ` is a root of `X²+X+1`, i.e.
`\barθ ∈ 𝔽₄ \ 𝔽₂` (a primitive cube root of unity); in particular `\barθ ∉ 𝔽₂`, so `v₂(θ) = 0`
(θ is a unit) and `\barθ ≠ \bar m` for every integer `m` (whose reduction lies in `𝔽₂`).
Therefore `m − θ` is a unit: `v₂(m − θ) = 0` for all `m ∈ ℤ`. Summing the non-special roots
gives `Σ_{θ ≠ ρ_b} v₂(m−θ) = 0`, and the only surviving root-term is `v₂(m − ρ_b)`. ∎

The constant `v₂(ℓ_b)` may be negative (e.g. `ℓ_6 = −1/90`, `v₂(ℓ_6) = −1`), reflecting that
`I_b` has integer-valued but non-integer coefficients; this is exactly the per-`b` offset
`c_b = v₂(ℓ_b)` observed numerically.

> **Corollary A1.** For an integer `m`, `I_b(m) = 0 ⟺ m ∈ {0,1,…,R}` or `m = ρ_b`. Since
> `m = ρ_b` forces `ρ_b ∈ ℤ`, the two-row d=4 law for `b` is equivalent to `ρ_b ∉ ℤ_{≥b}`.

**Proof.** `I_b(m) = 0 ⟺ v₂(I_b(m)) = +∞`. In Lemma A the first two terms are finite for any
fixed integer `m`; the right side is `+∞` iff `v₂(m − ρ_b) = +∞` (i.e. `m = ρ_b`) or one of the
forced factors vanishes (`m ∈ {0,…,R}`). The forced roots are `< b`; for `m ≥ b` the only way to
vanish is `m = ρ_b`. ∎

### Why this is the right reframing

Prior cycles recorded that `v₂(I_b(m))` is *unbounded* (it reaches 11 at `(b,m) = (6,18)`, and
grows without bound), and concluded "no finite 2-adic truncation closes the law." Lemma A
explains the unboundedness completely and benignly: it is **not** a pathology of central-trinomial
congruences but the single algebraic root `ρ_b` pulling the valuation up whenever an integer `m`
2-adically approximates it — `v₂(m − ρ_b)` is as large as the agreement of `m` with `ρ_b`. The
"spikes" of `v₂(I_b)` sit exactly at the integers nearest `ρ_b` in the 2-adic metric. So the
obstruction is not diffuse; it is one number `ρ_b`, and the law is the one statement
`ρ_b ∉ ℤ_{≥b}` — confirming and conceptually underpinning the single-candidate certificate of
06-09, with an exact valuation formula rather than a Newton-polygon estimate.

## 2. Supporting structure found this session (not load-bearing, but orienting)

**(a) Oscillatory (Chebyshev) form.** Writing `A = (1+αu)(1+βu)` with `αβ = 1`, `α+β = 1+i`,
and `α = e^{w}` (so `cosh w = (1+i)/2`, `w = p+iq` with `cosh p cos q = sinh p sin q = 1/2`):
```
    I_b(m) = Σ_k C(m,k) C(m, b−k) · sinh((2k−b)p) · sin((2k−b)q).
```
The factor `sin((2k−b)q)` is what makes `Q_b` carry many spread-out real roots (e.g. for
`b = 27` the real roots beyond `b` are ≈ `27.8, 42.2, 71.0, 162.1`). "Is one of these an integer"
is therefore a question with a transcendence flavour — consistent with the depth of the wall.

**(b) Uniform low 2-adic digits.** `ρ_b ≡ 2 (mod 16)` for all tested `b` (so `v₂(ρ_b) = 1` and
the next two digits are pinned); the digit at `2^4` and above is erratic. Four pinned bits is far
short of the `≈ 2log₂ b` bits one would need to force the low-digit integer out of `[b, U_b]`,
which is exactly why a uniform digit argument does not close `(★)`.

**(c) Half-integer near-roots.** For `b ≡ 2,3`, `Q_b(b − 1/2)` is a *unit fraction* `±1/N_b`
(tiny but nonzero); for `4 ∣ b` the same point `m = b − 1/2` is an **exact** root (the
half-integer factor divided out in the `4∣b` proof). This exhibits the `b mod 4` dichotomy as a
single phenomenon — the half-integer `b−1/2` is a genuine root exactly when `4 ∣ b` — but it
concerns half-integer, not integer, roots and so does not bear on `(★)`.

## 3. The gap, stated precisely (unchanged target, sharper footing)

> **Open `(★)`, uniform.** For all `b ≡ 2,3 (mod 4)`: the unique 2-adic root `ρ_b` is not a
> rational integer `≥ b`. Equivalently (Theorem 2): `Q_b` has no linear factor over `ℚ`.

Lemma A shows there is **no further 2-adic leverage** to extract: the valuation is *exactly*
`v₂(ℓ_b) + v₂((m)_{R+1}) + v₂(m − ρ_b)`, so any 2-adic certificate is precisely a statement that
the integer formed by the low `≈ 2log₂ b` digits of `ρ_b` avoids `[b, U_b]` — the digits are
erratic (no uniform pattern past `mod 16`), matching the established facts that no single witness
prime and no finite covering set exist. The gap is genuinely the **no-linear-factor**
(equivalently, for the empirically irreducible `Q_b`, irreducibility) of the parameter-Gegenbauer
family `g_b(m) = C_b^{(−m)}(−(1+i)/2)` — parameter-versus-degree, which the session brief
correctly flagged as the wall and instructed not to attack directly.

What Lemma A **does** give, concretely:
- the cleanest valuation certificate per `b` (exact, not estimated): finite `v₂(I_b(m))` for any
  integer `m ≠ ρ_b`, computed from `v₂((m)_{R+1})` and the 2-adic distance to `ρ_b`;
- a clean conceptual statement that the unbounded valuation is a feature (one root), not an
  obstruction internal to the congruence machinery — closing off the "central-trinomial
  congruence" route as a *source of new leverage* (Route α): there is none beyond `ρ_b` itself.

## 4. Verification

- Monic mod-2 factorisation `Q̂_b ≡ X(X²+X+1)^{(d−1)/2}` with `X²+X+1` irreducible: **all**
  `b ≡ 2,3 (mod 4)`, `6 ≤ b ≤ 80`. (`lemmaA_verify.monic_mod2_ok`.)
- Lemma A (single constant offset `c_b = v₂(ℓ_b)`, no fluctuation over `m`): **all**
  `b ≡ 2,3 (mod 4)`, `6 ≤ b ≤ 80`, each over `m ∈ [b, b+120)` (and spot-checked to `m = b+400`,
  including `b = 26,34,50,62`). (`lemmaA_verify.lemmaA_offset`.)
- `ρ_b ≡ 2 (mod 16)`; `ρ_b mod 32 ∈ {2,18}` (non-constant): `6 ≤ b ≤ 59`.

Files: `~/projects/code/lemmaA_verify.py` (verification harness, built on
`dfour_tworow.py`); this session reuses the prior `dfour_tworow` engine without modification.

## 5. Honest verdict

Did **not** close `(★)` uniformly. Produced the load-bearing lemma the brief asked for — an
**exact** closed form for `v₂(I_b(m))` that isolates the entire obstruction to the single 2-adic
distance `v₂(m − ρ_b)`, proving rigorously (given the established monic mod-2 factorisation) that
no central-trinomial/2-adic congruence can yield leverage beyond the location of `ρ_b` itself.
The residual gap is precisely `ρ_b ∉ ℤ_{≥b}`, i.e. no linear factor of `Q_b` — the
parameter-vs-degree irreducibility wall, unchanged in target but now resting on an exact
valuation identity rather than a Newton-polygon estimate.
