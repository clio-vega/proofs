# FINDINGS — the 𝔽₂ form of `I_b(m)` via the nilpotent truncation `s² ≡ 0`

**Date:** 2026-06-08 (code session)
**Scripts:** `dfour_tworow.py` (shared, independently cross-checked), `job1_fp2_recurrence.py`,
`job1_mod4_reduction.py`, `job1_v2_digits.py`
**Covers:** CODE.md Job 1.2 (headline), with Job 1.1 and Job 1.3 in the companion
`FINDINGS-job1-b0mod4.md`.

## Setup (reconstructed from scratch, not trusted from memory)

For the two-row shape `λ = (2m−b, b)`,
```
G(u; m) = (1−u)·(1 + s u + u²)^m,   s = 1+i,
I_b(m)  = Im( [u^b] G(u; m) )  ∈ ℤ[m],
I_b(m)  = (m)_{R+1} · Q_b(m),    (m)_{R+1} = m(m−1)···(m−R),  R = ⌊(b−1)/2⌋,
```
`q_b` = the monic primitive integer irreducible part of `Q_b`. The two-row d=4 law ⟺
`q_b` has no rational root `m ≥ b`. The module `dfour_tworow.py` builds `I_b, Q_b, q_b`
and **cross-checks `I_b` against an independent `u`-expansion** of `(1−u)(1+su+u²)^m` at
several integer `m` (passes for `b ≤ 16`; the `Q_b·(m)_{R+1} = I_b` reconstruction passes too).

## THE ALGEBRAIC REASON — nilpotent binomial truncation

Work in `R = ℤ[i]/(2)`. The key fact:
```
s = 1+i  ⟹  s² = (1+i)² = 2i ≡ 0  (mod 2).
```
**`s` is nilpotent mod 2.** Hence in `R[u]` the `m`-th power *truncates after two terms*
(every `(su)^n` with `n ≥ 2` contains `s² ≡ 0`):
```
(1 + s u + u²)^m  =  ((1+u²) + s u)^m
                  ≡  (1+u²)^m  +  m · s · u · (1+u²)^{m−1}      (mod 2).
```
Since `i ≡ s+1 (mod 2)`, a Gaussian integer `a+bi ≡ (a+b) + b·s`, so **`Im(·) =` the
`s`-coefficient.** The `s⁰` part `(1+u²)^m` is real and drops out of `Im`, leaving
```
I_b(m) ≡ [u^b]( (1−u)·m u (1+u²)^{m−1} )
       ≡ m·( [u^{b−1}] − [u^{b−2}] )(1+u²)^{m−1}      (mod 2).
```
Exactly one of `b−1, b−2` is even; the odd one contributes `0`, the even one `2j` gives
`C(m−1, j)`. With `R = ⌊(b−1)/2⌋` this collapses to a **closed form**:

> **(#)   `I_b(m) ≡ m · C(m−1, R)  (mod 2)`,  for every integer `m`, uniformly in `b`.**

This is a *proof skeleton*, not a finite check: the truncation argument is `b`-independent.

## What was verified (exact arithmetic)

- **(A) Truncation identity.** For `b ≤ 24` and several integer `m`, the reduced
  `[u^b]` Gaussian-integer coefficient of `(1−u)(1+su+u²)^m` equals that of the two-term
  truncation `(1+u²)^m + m s u (1+u²)^{m−1}`, mod 2 in both real and imaginary parts. ✓
- **(B) Closed form (#).** `I_b(m) − m·C(m−1,R)` has all **Mahler (forward-difference)
  coefficients ≡ 0 (mod 2)** for `b = 1…40`. (Mahler coefficients are the rigorous,
  finite certificate that two integer-valued polynomials agree mod 2 at *every* integer,
  not a sample of values.) ✓
- **(C) Downstream factorisation.** In `𝔽₂[m]`, for `b = 5…40`:
  ```
  q_b ≡ (m²+m+1)^{⌊(b−1)/4⌋}                 for b ≡ 0,1 (mod 4),
  q_b ≡ m · (m²+m+1)^{⌊(b−1)/4⌋}             for b ≡ 2,3 (mod 4).
  ```
  The `(m²+m+1)` multiplicity equals `⌊(b−1)/4⌋` throughout. ✓

## The decisive consequence for the proof

The extra factor `m` in the `b ≡ 2,3 (mod 4)` case is **exactly** why the 2-adic
certificate dies there: `q_b` has the root `m ≡ 0 (mod 2)`, so "no root mod 2" fails and
the 2-adic Newton/valuation argument cannot exclude an integer root. For `b ≡ 0,1 (mod 4)`
there is no such factor — `(m²+m+1)` is the irreducible `𝔽₂`-quadratic with **no root in
`𝔽₂`** — so `q_b` has no root mod 2, and the 2-adic argument closes those classes (already
proved for `b ≡ 1`, reduced to a mod-4 lemma for `b ≡ 0`; see companion findings).

So the `𝔽₂` trichotomy is now *explained*, not merely tabulated: it is the parity of `b`
deciding whether the surviving linear-in-`s` term `m·u·(1+u²)^{m−1}` lands its `u`-power on
an even index (→ `C(m−1,R)`, root-free quadratic part) or drags along the bare factor `m`.

## Honest scope

`(#)` and the truncation are genuinely uniform in `b` (proved, not sampled). The factorisation
`(C)` is verified `b ≤ 40`; its uniform form follows from `(#)` plus the standard `𝔽₂`
factorisation of `m·C(m−1,R)` (the `(m²+m+1)` powers come from `C(m−1,R) mod 2` via Lucas),
which the prove phase can finish. Nothing here touches the **hard** `b ≡ 2,3` non-vanishing
— that needs the odd-prime input of Job 2.
