# Two-row d=4 law, classes b ≡ 2,3 (mod 4): the unique 2-adic root, and a single-candidate certificate

**Date:** 2026-06-09 (prove session)
**Status:** (★) **not closed uniformly.** New (and now **proved unconditionally**, not merely
verified): the bare-factor-`m` "obstruction" makes `Q_b` have a **unique root `ρ_b` in ℤ₂**, so
`Q_b` has **at most one rational root**. This yields the cleanest known certificate (one candidate,
one prime, one evaluation per `b`), with which (★) is verified for all `b ≡ 2,3 (mod 4)` in
`[6, 171]`. The residual gap collapses to a single 2-adic statement: `ρ_b ∉ ℤ_{≥b}`.

---

## 0. The object and the gap

For `λ = (2m−b, b)`, `s = 1+i`, `A = 1+su+u²`, set `g_b(m)=[u^b]A^m`,
`G_b(m)=[u^b]((1−u)A^m)=g_b−g_{b−1}`, and `I_b(m) = Im G_b(m) ∈ ℤ[m]` (degree `b`).
The two-row d=4 fiber law `G_{(2m−b,b)}(i)=0 ⟺ (2,2)` reduces, after dividing out the forced
integer roots `m = 0,1,…,R` (`R=⌊(b−1)/2⌋`), to the monic-over-ℚ polynomial
`Q_b(m)` of degree `d = ⌊b/2⌋`, and to

> **(★)** for `b ≡ 2,3 (mod 4)`, `b ≥ 3`: `Q_b(m)` has no integer root `m ≥ b`.

Classes `b ≡ 0,1 (mod 4)` are already proved unconditionally (06-07/06-08, 2-adic valuation
identities). Closing (★) finishes the entire two-row d=4 law.

## 1. The engine (verified fresh against direct ℤ[i] coefficient extraction)

Writing `A = (su) + (1+u²)` and expanding,
`g_b(m) = Σ_{j} C(m,j) s^j C(m−j, (b−j)/2)` (only `j ≡ b mod 2`). Subtracting `g_{b−1}`,

> **Engine.** `I_b(m) = Σ_{j=1}^{b} (−1)^{b−j} Im(s^j) · C(m,j) · C(m−j, ⌊(b−j)/2⌋)`,
> with `Im(s^j) = 0` if `4 ∣ j`, and `Im(s^j) = {1,2,2}[j mod 4]·(−4)^{⌊j/4⌋}` otherwise.

Verified equal to the direct Gaussian-integer extraction `Im [u^b]((1−u)(1+su+u²)^m)` for all
tested `(b,m)`.

## 2. The mod-2 collapse, and why the 2-adic root is UNIQUE

> **Lemma 1 (mod-2 collapse).** `I_b(m) ≡ m · C(m−1, R) (mod 2)`, `R = ⌊(b−1)/2⌋`.

*Proof.* In the engine, `Im(s^j) = {1,2,2}[j mod 4]·(−4)^{⌊j/4⌋}`. For `j ≥ 4` the factor
`(−4)^{⌊j/4⌋}` is `≡ 0 (mod 2)`; for `j ∈ {1,2,3}`, `Im(s^j) = 1,2,2`, which mod 2 is `1,0,0`.
So mod 2 only the `j=1` term survives: `I_b ≡ (−1)^{b−1}·1·C(m,1)·C(m−1,R) ≡ m·C(m−1,R)`. ∎

The bare factor `m` is exactly what blocked every prior single-prime attempt: it produces the
root `m ≡ 0 (mod 2)`, so no congruence excludes an integer root. **But the same factor makes the
2-adic root unique.** The monic, 2-integral `Q_b` factors mod 2 as

> **Lemma 2 (mod-2 factorisation; prior cycle, Mahler coefficients; re-verified here for `6≤b≤80`).**
> ```
>    Q_b ≡ m · (m² + m + 1)^{(d−1)/2}   (mod 2),     m²+m+1 irreducible over 𝔽₂,
> ```
> where `d = ⌊b/2⌋` is odd for `b ≡ 2,3`, so `(d−1)/2 = ⌊(b−1)/4⌋ ∈ ℤ`.

*(Earlier confusion resolved: this clean form is for the **monic** `Q_b`. Reducing the **primitive
integer** model mod 2 instead corrupts it at `b=26,34,…` via an even leading coefficient — an
artifact, not a real break.)*

> **Theorem 2 (unique 2-adic root).** For every `b ≡ 2,3 (mod 4)`, `Q_b` has **exactly one** root
> `ρ_b` in `ℤ₂`, and `ρ_b` is even (`v₂(ρ_b) ≥ 1`).

*Proof.* By Lemma 2, over `𝔽₂` the only root of `Q_b` is `m=0` (since `m²+m+1` has no `𝔽₂`-root:
its values at `0,1` are `1`), and it is **simple** — `(Q_b mod 2)'(0) = (m²+m+1)^{(d−1)/2}|_{m=0}
= 1 ≠ 0`. Hensel's lemma lifts this simple root to a **unique** `ρ_b ∈ ℤ₂` with `ρ_b ≡ 0 (mod 2)`.
Conversely any root of `Q_b` in `ℤ₂` reduces mod 2 to an `𝔽₂`-root, hence to `m=0`, hence (Hensel
uniqueness) equals `ρ_b`. So `ρ_b` is the unique `ℤ₂`-root, and `ρ_b ≡ 0 (mod 2)`. ∎

(The 2-adic Newton polygon `[slope −1, run 1] ⊕ [slope 0, run d−1]`, computed for `b ≤ 31`, sharpens
this to `v₂(ρ_b) = 1` and shows the `d−1` unit roots lie in the unramified quadratic extension — but
only `v₂(ρ_b) ≥ 1` is needed below.)

> **Corollary 3.** `Q_b` has **at most one rational root**, namely `ρ_b` (and it is even).

*Proof.* A rational root is an integer (monic), hence lies in `ℤ ⊂ ℤ₂` and is a root of `Q_b`
there; by Theorem 2 it must equal the unique 2-adic root `ρ_b`. ∎

This **inverts** the prior reading of the obstruction: the bare factor `m` is not a wall, it is the
reason the candidate set collapses from `~0.1 b²` integers to a **single** number `ρ_b`, computable
from the **single prime 2** to any precision.

## 3. The single-candidate certificate (rigorous, per b)

> **Certificate.** Fix `b ≡ 2,3 (mod 4)`. Let `U_b = 1 + max_i |a_i|` be the Cauchy bound on the
> roots of the monic `Q_b = Σ a_i m^i` (all real roots satisfy `|m| ≤ U_b`). Choose `K` with
> `2^K > U_b` and Hensel-lift the 2-adic root to `r := ρ_b mod 2^K ∈ [0, 2^K)`.
> Then the **only** possible integer root `m_0 ∈ [b, U_b]` is `m_0 = r`. If `r ∉ [b, U_b]`, or
> `Q_b(r) ≠ 0`, then `Q_b` has **no** integer root `≥ b`, i.e. (★) holds for `b`.

*Proof.* An integer root `m_0 ≥ b` satisfies `0 < m_0 ≤ U_b < 2^K`, and `m_0 = ρ_b` in `ℤ₂`
(Corollary 3), so `m_0 ≡ ρ_b (mod 2^K)`; the unique representative in `[0,2^K)` is `r`, whence
`m_0 = r`. So `r` is the sole candidate; one exact evaluation `Q_b(r)` settles it. ∎

This needs **one** prime, **one** lifted root, **one** integer evaluation per `b` — far cheaper than
rational-root enumeration / full factorisation of a degree-`d` polynomial.

> **Theorem 4 (frontier).** (★) holds for **all** `b ≡ 2,3 (mod 4)` with `6 ≤ b ≤ 171`.
> Consequently the two-row d=4 law holds for these `b`.

Cross-checked against exact rational-root enumeration (`ground_roots`) for `b ≤ 40`: agreement.
(Examples: `b=6`, `U_6=54`, sole candidate `r=18`, `Q_6(18)=1920≠0`; `b=7`, `U_7=791`, sole
candidate `r=610`, `Q_7(610)≠0`. In both the single evaluation excludes the lone candidate.)

## 4. The residual gap, stated precisely

> **Open ((★), uniform).** For all `b ≡ 2,3 (mod 4)`: the unique 2-adic root `ρ_b ∈ ℤ₂` is **not**
> a rational integer in `[b, ∞)`. Equivalently: the low `≈ 2 log₂ b` 2-adic digits of `ρ_b`, read
> as an integer, do not land in `[b, U_b]` (or, if they do, `Q_b` does not vanish there).

By Theorem 2 / Corollary 3 this is now a statement about **one** 2-adic number per `b`, not about a
polynomial's whole integer-root locus. Two equivalent shapes of the same gap:

- (digit form) `ρ_b`'s 2-adic expansion is not eventually `0` and not eventually `1` within its
  first `≈2 log₂ b` digits — i.e. `ρ_b ∉ ℤ` in the relevant window;
- (Galois form) `Q_b` has no linear factor over `ℚ`. For irreducible `Q_b` (empirically all of
  them, `b ≤ 171`) this is automatic; by **Jordan's derangement theorem + Chebotarev** an
  irreducible `Q_b` of degree `≥ 2` has a positive density of primes `p` with `Q_b` rootless mod
  `p`, each a finite certificate. So the entire uniform gap = **irreducibility (or just
  no-linear-factor) of the Gegenbauer-in-parameter family** `g_b(m) = C_b^{(−m)}(−(1+i)/2)`.

What is now **rigorously dead** (do not re-attempt): single-prime Newton polygon (all slopes
integral at `p=2`; never a full non-integral run at `p ≤ 7`), uniform witness prime (irregular:
`7,11,11,17,17,29,19,41,…`), finite covering set, analytic/log-concavity domination (real roots
are in-range, 06-08 Prop 3), `q`-Lucas / Gaussian-binomial shortcut (`G_b` is not that object).

## 5. Files

`~/projects/scratch/`:
`qtest.py, factor.py, padic.py, witness.py, newton.py, newton2.py` (route-death, fresh);
`uniq2adic.py, uniq2.py, mod2fact.py, rho.py` (unique 2-adic root, Theorem 2);
`cert_fast.py, cert_push.py` (the certificate, Theorem 4); notebook
`prove-20260609-tworow-d4-b23.md`.

## Summary

Not a closure. But a genuine reframing: the bare-factor-`m` that defeated the local-prime program
is exactly what makes the 2-adic root **unique**, collapsing (★) to a **single-candidate** test
that is rigorous per `b` and verifies the law to `b ≤ 171` from one prime. The uniform gap is now a
single sharp statement — `ρ_b ∉ ℤ_{≥b}`, i.e. no linear factor of the parameter-Gegenbauer family.
