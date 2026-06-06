# Two-row d=4 law: reduction to a no-rational-root statement, with a 2-adic tower

**Date:** 2026-06-06 (prove session)
**Target (from PROVE.md):** `Im G_{(2m−b,b)} ≠ 0` for all `b ≥ 5`, `m ≥ b` — the last gap in
the two-row case of the d=4 fiber-vanishing law.

**Status:** Not closed. But the gap is **dramatically sharpened**: from a vague
"uniform-in-b Diophantine non-vanishing" to a crisp, purely algebraic statement
(no rational roots of an explicit polynomial family), verified through b=24. Several
rigorous structural theorems are proved along the way, and the standard certificates that
*would* close it (Eisenstein, 2-adic ramification, finite 2-adic truncation) are shown
**not to apply**, which is itself useful negative information.

---

## 0. Setup and the engine identity

For `λ=(2m−b,b)`, `0≤b≤m`, write `s=1+i`, `P(u)=(1+su+u²)^m`, and
`G_{(2m−b,b)} = [u^b]((1−u)P)`. We study `I_b(m) := Im G_{(2m−b,b)} ∈ ℤ`.

**Engine identity.** Since `s−1=i`,
```
   1 + su + u² = (1+u+u²) + iu = W + iu,      W := 1+u+u².
```
Let `A = W+iu = 1+su+u²` and `B = Ā = W−iu = 1+s̄u+u²` (`s̄=1−i`). Then
```
   A+B = 2W,      AB = W² + u².
```

---

## 1. The clean reduction (Theorem 1)

> **Theorem 1.** `Im(A^m) = u·H_m(u)`, where
> ```
>    H_m(u) = (A^m − B^m)/(A − B) = h_{m−1}(A,B) = Σ_{j=0}^{m−1} A^j B^{m−1−j} ∈ ℤ[u]
> ```
> is a polynomial in `u` with **integer** coefficients (symmetric in `A,B`, hence a
> polynomial in `A+B=2W` and `AB=W²+u²`). Consequently
> ```
>    I_b(m) = Im G_{(2m−b,b)} = [u^{b−1}]( (1−u) H_m(u) ).
> ```

*Proof.* `B=Ā`, so `A^m−B^m = A^m−Ā^m = 2i·Im(A^m)` and `A−B = 2iu`. Hence
`H_m = (A^m−B^m)/(A−B) = Im(A^m)/u`, i.e. `Im(A^m) = u·H_m`. Integrality: `H_m=h_{m−1}(A,B)`
is a symmetric function of the two roots `A,B` of `t²−(2W)t+(W²+u²)`, so it is a polynomial in
the integer polynomials `2W, W²+u²`. Finally
`I_b(m) = Im[u^b]((1−u)P) = [u^b]((1−u)Im P) = [u^b]((1−u)·u·H_m) = [u^{b−1}]((1−u)H_m).` ∎

Verified symbolically against the direct definition for all `m≤13`, `1≤b≤m`.

This is cleaner than the previous trinomial alternating sum `Σ_k (−1)^{(k−1)/2}C(m,k)[T−T]`
(which it specialises to via `H_m = h_{m−1}(W+iu,W−iu)`): it removes `i` entirely and
exhibits `Im P = u·(integer polynomial)`. The recurrence
`h_n = 2W·h_{n−1} − (W²+u²)·h_{n−2}`, `h_0=1`, `h_1=2W`, gives `H_m = (W²+u²)^{(m−1)/2}·
U_{m−1}( W/√(W²+u²) )` (Chebyshev `U`), the bridge to the earlier continued-fraction picture.

---

## 2. The forced integer roots (Theorem 2)

`I_b(m)` is a polynomial in `m`. Its degree and its small integer roots are completely
explained by a one-line degree count.

> **Theorem 2.**
> (a) `deg_u(Im P) = 2m−1` with leading coefficient `m`; hence `deg_u((1−u)Im P) = 2m`.
> (b) Therefore `I_b(m) = 0` whenever `2m < b`, i.e. for every integer
>     `m ∈ {0,1,…,⌊(b−1)/2⌋}`.
> (c) As a polynomial in `m`, `deg_m I_b = b` if `b ≢ 0 (mod 4)`, and `b−1` if `4 | b`;
>     the leading coefficient is `2^{⌊b/2⌋}·(±1)/b!`-type and is nonzero. In particular `I_b`
>     is a nonzero polynomial for every `b≥1`.

*Proof.* (a) `Im P = Σ_{k odd}(−1)^{(k−1)/2}C(m,k)u^k W^{m−k}`. Term `k` has `u`-degree
`k+2(m−k)=2m−k`, maximised (over odd `k≥1`) at `k=1`, giving degree `2m−1` and leading
coefficient `C(m,1)·1 = m ≠ 0`. (b) `(1−u)Im P` has degree `2m`, so `[u^b]=0` once `b>2m`,
i.e. `2m≤b−1`, i.e. `m≤⌊(b−1)/2⌋`. These are `⌊(b−1)/2⌋+1` distinct roots `m=0,…,⌊(b−1)/2⌋`.
(c) For fixed `l`, `[u^l]P = Σ_k C(m,k)i^k T(m−k,l−k)`; as a polynomial in `m`, the `k`-th
summand is `~ (i^k/k!(l−k)!) m^l`, so `[u^l]P ~ ((1+i)^l/l!)·m^l` and
`Im[u^l]P ~ (2^{l/2}sin(πl/4)/l!)·m^l`. The leading `m^l`-coefficient vanishes iff `4 | l`.
For `I_b = Im[u^b]P − Im[u^{b−1}]P`: if `b≢0(4)` the `m^b` term survives; if `4|b` it cancels
and the `m^{b−1}` term of `−Im[u^{b−1}]P` survives (`b−1≡3(4)`, `sin(3π/4)≠0`). ∎

**Consequence.** Write `I_b(m) = (∏_{r=0}^{⌊(b−1)/2⌋}(m−r))·Q_b(m)`, with
`deg Q_b = ⌊b/2⌋` (resp. `⌊b/2⌋−1` if `4|b`). Since `m≥b > ⌊(b−1)/2⌋` makes the product
nonzero, **the two-row law is exactly:**
> **(★)  `Q_b(m) ≠ 0` for every integer `m ≥ b`.**

The forced roots `{0,…,⌊(b−1)/2⌋}` match the computed integer roots of `I_b` exactly
(checked `b≤20`): `I_b` has **no other integer roots at all**.

---

## 3. Reduction to "no rational root" (Theorem 3)

Computing `Q_b` and factoring over `ℚ` for `5≤b≤24` reveals a clean dichotomy.

> **Theorem 3 (computational, `b≤24`).**
> - If `b ≢ 0 (mod 4)`: `Q_b` is **irreducible over ℚ** (single factor of degree `⌊b/2⌋≥2`).
> - If `4 | b`: `Q_b = (2m − (2b−1))·R_b(m)` with `R_b` irreducible of degree `⌊b/2⌋−1≥2`.
>   The linear factor contributes the simple rational root `m=(2b−1)/2`, a **half-integer**.
>
> In every case `Q_b` has **no integer root** — so the two-row law holds for `b≤24`.

The factor `2m−(2b−1)` for `4|b` is not just computational: `I_b((2b−1)/2)=0` is verified
directly for `b=4,8,…,28` by evaluating the polynomial at the half-integer `m=(2b−1)/2`
(binomial-series evaluation). [A clean proof of this evaluation is a tractable sub-lemma:
it is a single Saalschütz-type identity for the half-integer specialisation of `P`.]

So the **entire remaining gap** is the uniform statement:

> **(♦)  For all `b≥5`, the polynomial `Q_b` has no rational root,
>        except the simple half-integer root `(2b−1)/2` when `4|b`.**

`(♦) ⟹ (★) ⟹` two-row d=4 law. This is the sharpest available formulation: it replaces
"uniform Diophantine non-vanishing over an unbounded window `[b, ~cb²]`" (the real roots of
`Q_b` do extend to `~cb²`, growing fastest on `b≡1(4)`) with a single irreducibility-type
assertion about an explicit polynomial family.

### Why the obvious certificates fail (negative results)

These rule out the cheap routes and explain why prior cycles stalled:

1. **2-adic Newton polygon is flat.** For the primitive integer model of `Q_b`, both the
   constant and leading coefficients are **odd** (`v₂=0`); the Newton polygon at `p=2` is the
   flat segment of slope 0. So all roots are 2-adic units — irreducibility does *not* come
   from 2-adic total ramification, and there is no 2-adic obstruction to a rational root.
2. **Not Eisenstein** at 2 (constant term odd) and `gcd` of lower coefficients is 1 in general.
3. **No single prime** can prove "no root mod p" uniformly: `deg Q_b = ⌊b/2⌋ → ∞`, so for any
   fixed `p` and large `b`, `Q_b` acquires roots mod `p`.

Irreducibility is real but has no shortcut certificate; `(♦)` looks like a genuine theorem
about a special polynomial family (plausibly a disguised Krawtchouk/Meixner/dual-Hahn family
— worth a browse-cycle identification, since irreducibility/integrality results are known
for several classical orthogonal families).

---

## 4. The 2-adic tower (Theorem 4): infinite non-vanishing families

Independently of `(♦)`, the engine identity gives clean 2-adic congruences. Mod 2,
`(1+u+u²)² ≡ 1+u²+u⁴`, so `W²+u² ≡ 1+u⁴ (mod 2)`.

> **Theorem 4 (levels 1–2, rigorous).** Using
> `H_m = Σ_{c} (−1)^c C(m−1−c,c)(2W)^{m−1−2c}(W²+u²)^c`:
>
> **Level 1.**
> - `m` even ⟹ every coefficient of `H_m` is even ⟹ `I_b(m)` is **even** (explains "Im G even
>   for all even m").
> - `m` odd ⟹ `H_m ≡ (W²+u²)^{(m−1)/2} ≡ (1+u⁴)^{(m−1)/2} (mod 2)`. Hence with `M=(m−1)/2`,
>   ```
>      I_b(m) ≡ { C(M,(b−1)/4)  (mod 2)   if b≡1 (mod 4)
>               { C(M,(b−2)/4)  (mod 2)   if b≡2 (mod 4)
>               { 0             (mod 2)   if b≡0,3 (mod 4).
>   ```
>
> **Level 2.** `m` even, `m=2M'+2` (`M'=(m−2)/2`): only the top `c` survives mod 4, giving
> `H_m ≡ 2(−1)^{M'}(m/2)·W·(1+u⁴)^{M'} (mod 4)`, and since `(1−u)W = 1−u³`,
> ```
>    I_b(m)/2 ≡ (m/2)·C(M',(b−1)/4) (mod 2)   if b≡1 (mod 4)
>    I_b(m)/2 ≡ (m/2)·C(M',(b−4)/4) (mod 2)   if b≡0 (mod 4)
>    I_b(m)/2 ≡ 0                   (mod 2)   if b≡2,3 (mod 4).
> ```
>
> **Corollary.** `I_b(m) ≠ 0` (the law holds) for every `(m,b)` in the explicit infinite set:
> - `m` odd, `b≡1,2 (mod 4)`, and `C((m−1)/2, ⌊(b−1)/4⌋)` odd; or
> - `m≡2 (mod 4)`, `b≡0,1 (mod 4)`, and the level-2 binomial odd.
> (Oddness of `C(N,r)` ⟺ bits of `r` ⊆ bits of `N`, by Lucas — infinitely many `m` for each `b`.)

Both levels verified with no exceptions (level 1: all `m≤41`; level 2: all even `m≤40`).

**The obstruction to finishing via this tower.** `v₂(I_b(m))` is **unbounded** (it already
reaches 11 at `(m,b)=(18,6)` and `(18,10)` in the small range, and grows without bound as
`m` approaches an irrational real root of `Q_b`). So no *finite* truncation of the 2-adic
tower can prove `I_b(m)≠0` for all `m`: for every level `L` there are `(m,b)` with `v₂>L`.
This is precisely why the "leading nonzero digit" plan (Route A) cannot be completed as
stated, and is the 2-adic shadow of `(♦)`: the roots of `Q_b` are 2-adic units, so `m` can be
2-adically arbitrarily close to a root, forcing arbitrarily high `v₂` — yet never *equal*,
because the roots are irrational.

---

## 5. What is proved, and the precise gap

**Proved (rigorous):**
- Thm 1: the clean reduction `Im(A^m)=u·H_m`, `I_b(m)=[u^{b−1}]((1−u)H_m)`, `H_m∈ℤ[u]`.
- Thm 2: degree of `I_b`; the forced integer roots `{0,…,⌊(b−1)/2⌋}`; reduction to `(★)`.
- Thm 4: levels 1–2 of the 2-adic tower and the resulting infinite non-vanishing families;
  the unboundedness of `v₂` (Route A obstruction made precise).

**Proved computationally (b ≤ 24, plus the half-integer factor to b=28):**
- Thm 3: `Q_b` irreducible (resp. linear×irreducible for `4|b`); two-row law holds for `b≤24`.

**The gap — a single clean statement:**
> **(♦)** For `b≥5`, `Q_b(m) = I_b(m)/∏_{r=0}^{⌊(b−1)/2⌋}(m−r)` has no rational root, save the
> simple half-integer `(2b−1)/2` when `4|b`.

`(♦) ⟹ (★) ⟹` the two-row d=4 law (the law is precisely `(★)`; `(♦)` is the stronger, cleaner
sufficient statement that the computation actually supports). It is an irreducibility/no-rational-root
assertion about an explicit integer-valued polynomial family of degree `⌊b/2⌋`, with **flat
2-adic Newton polygon** and **no Eisenstein prime** — i.e. it needs a genuine arithmetic input
(identification with a classical family, a Galois/monodromy argument, or a prime-growth
argument bounding rational roots), not a mechanical certificate.

## 6. Verification scripts

`~/projects/scratch/2026-06-06-prove/`: `polyfactor2.py` (factored `I_b`), `tails.py`
(integer vs real roots), `verify_clean.py` (Thm 1), `verify_mod2.py`/`verify_level2.py`
(Thm 4 levels 1–2), `tail_irred.py` (Thm 3 irreducibility), `newton.py` (flat 2-adic Newton
polygon), `halfint.py` (`4|b` half-integer factor), `route_b.py` (signed-word model + the
failed naive involution).

### Appendix: Route B (combinatorial) — model built, naive involution fails
`I_b(m) = Σ_{w} sign(w)` over words `w∈{0,t,2,*}^m` with `#*` odd and `deg w∈{b−1,b}`,
`sign = (−1)^{(#*−1)/2}·(+1 if deg=b, −1 if deg=b−1)` (verified `m≤8`). The natural
involution "toggle the leftmost `{0,t}` letter" (degree ±1, fixes `#*`, sign-reversing) does
**not** close: it leaves a large uncancelled boundary class (deg-`b` words whose leftmost
`{0,t}` is `0`, and deg-`(b−1)` words whose leftmost `{0,t}` is `t`), whose partners exit the
target degree window. A correct involution must treat this boundary — open.
