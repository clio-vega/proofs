# FINDINGS — Job 2: exact central-trinomial identity for I_b(m)

**Date:** 2026-06-09 (code session)
**Script:** `job2_trinomial.py`
**Verified:** exact polynomial identity in ℤ[m] for all b ≤ 50.

## The clean factorisation that drives everything

With s = 1 + i,

  **1 + s u + u² = (1 + u + u²) + i u = W + i u**,   W := 1 + u + u².

(Verified: difference is identically 0.) Therefore

  (1 + s u + u²)^m = (W + i u)^m = Σ_k C(m,k) iᵏ uᵏ W^{m−k},

and taking the imaginary part keeps only **odd k** (iᵏ = i·(−1)^{(k−1)/2}):

  Im (1+su+u²)^m = Σ_{k odd} (−1)^{(k−1)/2} C(m,k) uᵏ W^{m−k}.

## The exact index map (TRI)

Applying [u^b]( (1−u) · ) to the above gives the explicit signed combination
the PROVE session needs:

> **(TRI)**  I_b(m) = Σ_{k odd, k ≤ b} (−1)^{(k−1)/2} C(m,k)
>                       · [ τ(m−k, b−k) − τ(m−k, b−1−k) ],
>
> where τ(n, j) = [u^j](1+u+u²)ⁿ = Σ_l C(n,l) C(n−l, j−2l) is the trinomial
> coefficient (triangle **A027907**).

**Verified exactly as polynomials in ℤ[m] for every b ≤ 50** (`(TRI) holds: True`).

## Where A002426 enters — the central diagonal

The trinomial coefficient τ(n, j) is **central** exactly when j = n, and the
central trinomial coefficients are

  **T_n := τ(n, n) = A002426 = 1, 1, 3, 7, 19, 51, 141, 393, 1107, 3139, …**

GF 1/√(1 − 2x − 3x²), recurrence (n+1)T_{n+1} = (2n+1)T_n + 3n T_{n−1} (both
re-verified). In (TRI), the argument pair (m−k, b−k) is central precisely on the
**diagonal m = b**: there the k=1 term contributes C(b,1)·[T_{b−1} − τ(b−1,b−2)],
i.e. the central trinomial T_{b−1} is the leading trinomial coefficient. So
"Im G = (−1)ⁿ A002426" is the **diagonal m = b specialisation** of (TRI), with the
full off-diagonal object being the trinomial *triangle* A027907, not a single
A-sequence.

## The m = b diagonal sequence

I_b(b) for b = 1..15:  1, 0, 2, 4, 6, 16, 20, 8, −98, −608, −2388, −8040,
−24244, −66208, −162552.  (No clean OEIS match was confirmable offline; the
sign change at b=9 reflects that this is *not* purely the positive central
trinomial but the signed (1−u)-difference combination above. Recorded for the
record; not load-bearing.)

## Use to PROVE

(TRI) is the exact, symbolic-in-m expression of I_b(m) via trinomial
coefficients. It is the cleanest available closed form and the natural place to
read off arithmetic (e.g. the W = 1+u+u² structure makes the mod-2 reduction
q_b ≡ m·(m²+m+1)ᵏ transparent: W ≡ 1+u+u² and the central-trinomial recurrence
mod 2 is governed by m²+m+1). The D-finite/holonomic nature of τ(·,·) (P-finite
in both indices) means I_b(m) is holonomic in m for each fixed b.
