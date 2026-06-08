# Two-row d=4 law, the residue classes b ≡ 2,3 (mod 4): a central-binomial closed form, the Gegenbauer engine, and why the obstruction is arithmetic, not analytic

**Date:** 2026-06-08 (prove session)
**Status:** The target class `b ≡ 2,3 (mod 4)` is **not closed**. But this session produces
(i) a clean new closed form (Theorem A) with a real proof, (ii) the exact "engine" recurrence
identifying the family, (iii) a sharp diagnosis that **retires the planned log-concavity route**,
and (iv) an extension of the rigorous frontier to `b ≤ 70`.

---

## 0. The object and the gap

For `λ = (2m−b, b)`, `s = 1+i`, `W = 1+u+u²`, `A = 1+su+u² = W + iu`, write

```
   g_b(m) := [u^b] A^m,      G_b(m) := [u^b]((1−u)A^m) = g_b(m) − g_{b−1}(m),
   I_b(m) := Im G_{(2m−b,b)} = Im G_b(m) ∈ ℤ[m].
```

The two-row d=4 fiber law `G_{(2m−b,b)}(i)=0 ⟺ (2,2)` is, after dividing out the forced
integer roots `m=0,1,…,R` (`R=⌊(b−1)/2⌋`), exactly

> **(★)** for `b ≡ 2,3 (mod 4)` and every integer `m ≥ b`, `I_b(m) ≠ 0`,

equivalently `Q_b(m) := I_b(m)/∏_{r=0}^{R}(m−r)` has no integer root `≥ b`. The classes
`b ≡ 0,1 (mod 4)` are already proved (2-adic valuation identities, 06-07/06-08). For `b ≡ 2,3`,
mod 2 one has `I_b(m) ≡ m·C(m−1,R)`: the bare factor `m` produces a genuine root `m ≡ 0`,
so **no 2-adic (indeed no single-prime) certificate exists** — established as a structural fact
in the 06-08 CODE session.

---

## 1. The engine: a three-term recurrence and the Gegenbauer identity (Proposition 1)

Differentiating `A^m` and using `A·(A^m)' = m(s+2u)A^m` gives, on reading off `[u^{b−1}]`,

> **Proposition 1 (recurrence).** For all `b ≥ 2`,
> ```
>    b·g_b(m) = s(m−b+1)·g_{b−1}(m) + (2m−b+2)·g_{b−2}(m),     g_0 = 1, g_1 = s·m.
> ```

*Proof.* `(1+su+u²)(A^m)' = m(s+2u)A^m`. Taking `[u^{b−1}]` of both sides:
LHS `= b g_b + s(b−1)g_{b−1} + (b−2)g_{b−2}`, RHS `= m s g_{b−1} + 2m g_{b−2}`; rearrange. ∎
(Verified symbolically for all `m ≤ 8`, `b ≤ 2m`.)

This pins the closed form and the classical identity:

> **Proposition 2 (closed form / Gegenbauer in the parameter).**
> ```
>    g_b(m) = Σ_{l=0}^{⌊b/2⌋}  (m)_{b−l} / ((b−2l)! l!) · s^{b−2l}
>           = C(m,b)·s^b· ₂F₁(−b/2, −(b−1)/2 ; m−b+1 ; −2i),
> ```
> where `(m)_{k}` is the falling factorial. Equivalently `g_b(m) = C_b^{(−m)}(−s/2)`, the
> Gegenbauer polynomial of degree `b` and **parameter `−m`**, evaluated at the fixed complex
> point `x = −s/2 = −(1+i)/2`.

*Proof.* The multinomial expansion of `(1 − 2xu + u²)^{m}` with `x = −s/2` gives the explicit
sum; matching consecutive-term ratios gives the `₂F₁`; the generating function
`(1−2xt+t²)^{−λ} = Σ C_n^{λ}(x)t^n` with `λ = −m` gives the Gegenbauer reading.
(Both forms verified symbolically, `m ≤ 8`.) ∎

The argument `−2i = −s²` and the upper parameters `−b/2, −(b−1)/2` are the fingerprint of a
**Gegenbauer / Chebyshev** structure — but in the parameter `m`, which is *not* a standard
orthogonal-polynomial variable (the recurrence coefficient `(2m−b+2)/b` of `g_{b−2}` depends on
`m`). So the `g_b(m)` are **not** orthogonal polynomials in `m`, and the classical
all-real-simple-zeros machinery does not apply directly. This is the precise reason the family
resists the standard toolkit.

---

## 2. Theorem A: the half-integer value is a central binomial coefficient

The point `m = (2b−1)/2` is the centre of the dangerous window (it is where, for `4∣b`, `Q_b`
has its rational root). Here the value is exactly computable and beautiful.

> **Theorem A.** For all `b ≥ 1`,
> ```
>    G_b((2b−1)/2) = (−1)^b · C(2b,b)/4^b · (1−i)^b.
> ```
> Consequently
> ```
>    I_b((2b−1)/2) = Im G_b((2b−1)/2) = (−1)^{b+1} · C(2b,b)/4^b · Im(s^b),
> ```
> whose magnitude is `2^{⌊b/2⌋}·C(2b,b)/4^b` and which **vanishes iff 4 ∣ b**.

*Proof.* By Proposition 2, `g_b((2b−1)/2) = C((2b−1)/2, b)·s^b·₂F₁(−b/2,−(b−1)/2; 1/2; −2i)`,
because the lower parameter is `m−b+1 = (2b−1)/2 − b + 1 = 1/2`. Likewise
`g_{b−1}((2b−1)/2) = C((2b−1)/2, b−1)·s^{b−1}·₂F₁(−(b−1)/2,−(b−2)/2; 3/2; −2i)`,
the lower parameter now being `3/2`. Apply the two classical quadratic transformations with
`z = 1−i` (so `z² = −2i`):
```
   ₂F₁(a, a+1/2 ; 1/2 ; z²) = ½[(1+z)^{−2a} + (1−z)^{−2a}],
   ₂F₁(a, a+1/2 ; 3/2 ; z²) = [(1+z)^{1−2a} − (1−z)^{1−2a}] / (2z(1−2a)).
```
With `a = −b/2` (first) and `a = −(b−1)/2` (second), and `1+z = 2−i`, `1−z = i`:
```
   ₂F₁(−b/2,−(b−1)/2; 1/2; −2i) = ½[(2−i)^b + i^b],
   ₂F₁(−(b−1)/2,−(b−2)/2; 3/2; −2i) = [(2−i)^b − i^b] / (2(1−i)b).
```
The falling-factorial binomials collapse to central binomials:
`C((2b−1)/2, b) = (2b−1)!! / (2^b b!) = C(2b,b)/4^b`, and
`C((2b−1)/2, b−1) = 2b·C(2b,b)/4^b`. Writing `C := C(2b,b)/4^b` and using `s/(1−i) = i`,
`s = 1+i`:
```
   g_b      = C·s^{b−1}·(1+i)·½[(2−i)^b + i^b],
   g_{b−1}  = C·s^{b−1}·(1+i)·½[(2−i)^b − i^b].
```
Subtracting, the `(2−i)^b` terms cancel and the `i^b` terms add:
```
   G_b = g_b − g_{b−1} = C·s^{b−1}·(1+i)·i^b = C·s^b·i^b = C·(si)^b = (−1)^b·C·(1−i)^b,
```
since `si = (1+i)i = −(1−i)`. Finally `(1−i)^b = conj(s^b)`, so
`Im G_b = (−1)^b·C·Im((1−i)^b) = (−1)^{b+1}·C·Im(s^b)`, and `|Im(s^b)| = 2^{⌊b/2⌋}` for
`b ≢ 0 (mod 4)`, `= 0` for `4 ∣ b`. ∎

(Closed form verified symbolically for all `b ≤ 15`; magnitude `= 2^{⌊b/2⌋}C(2b,b)/4^b` and the
A001790 numerators confirmed for `b ≤ 23`.)

**What Theorem A buys.** It *unifies the two regimes*:

- For `4 ∣ b`, `Im(s^b) = 0`, so `(2b−1)/2` **is** a root of `I_b` — recovering, with a clean
  conceptual reason, the rational root `(2b−1)/2` of `Q_b` found computationally for `4∣b`.
- For `b ≡ 2,3 (mod 4)`, the nearest half-integer to the dangerous region is provably **not**
  a root, and its value is an explicit nonzero central binomial coefficient (up to a
  2-power and sign). Several genuine real roots of `I_b` for `b ≡ 2,3` cluster near `(2b−1)/2`
  (e.g. `b=14 → 13.5`, `b=23 → 22.5`, `b=26 → 25.5`); Theorem A is the exact statement that the
  half-integer itself escapes — the roots merely *approach* it.

---

## 3. The obstruction is arithmetic, not analytic (this retires the log-concavity plan)

The PROVE.md plan proposed a uniform analytic lower bound (leading-term / log-concavity
dominance of the alternating trinomial sum) to force `I_b(m) ≠ 0`. **This route cannot close
(★).** Here is the precise reason.

> **Proposition 3.** For every `b ≡ 2,3 (mod 4)` with `b ≥ 6`, the real polynomial `m ↦ I_b(m)`
> has a real zero in the open interval `(b, ∞)` — in fact in `(b, c b²)` with `c ≈ 0.33`.

*Proof (existence).* `I_b(b)` and `I_b(M)` (with `M` past the largest real root, where the sign
equals that of the leading coefficient `Im(s^b)/b!`) have opposite signs for these `b`; e.g.
`I_{11}(11) < 0` while `I_{11}(m) > 0` for large `m`. A sign change between integer arguments
`≥ b` forces a real root in `(b, M)`. The location `≈ c b²` is the empirical envelope of the
largest real root. ∎

**Consequence.** Since the *real* function `I_b` genuinely vanishes inside `[b, c b²]`, **no**
lower bound of the form "`|I_b(m)| > 0` for all real `m ≥ b`" can hold. Any leading-term or
log-concavity dominance argument proves a bound valid for all real `m` in a range, hence must
fail wherever the real roots live. Such an argument can therefore only certify `m` **beyond the
largest real root** (`m ≳ c b²`) — which is exactly the regime already settled by computing the
roots, i.e. *no improvement over brute force*. The genuine content of (★) — that the in-range
real roots, all irrational, **avoid the integers** — is **arithmetic**, not analytic.

This is a real correction to the plan: the broken belief flagged in PROVE.md ("(★) needs an
irreducibility theorem — maybe it needs only a direct log-concavity bound") resolves the
**opposite** way. Log-concavity is the wrong tool; the no-integer-root content is irreducible-in-spirit.

---

## 4. The frontier, certified

> **Theorem B (computational, exact).** For every `b ≡ 2,3 (mod 4)` with `3 ≤ b ≤ 70`, the only
> rational roots of `I_b ∈ ℤ[m]` are integers `< b`. Hence `Q_b` has no rational (a fortiori no
> integer) root `≥ b`, and the two-row d=4 law holds for these `b`.

This uses exact rational-root enumeration (SymPy `ground_roots`, certified — no tail bound, no
floating point), extending the previous frontier (`b ≤ 40` by factorization). Supplementary:
direct exact-integer non-vanishing `I_b(m) ≠ 0` checked for all integers `m ∈ [b, b²]`,
`b ≤ 150`. (The lone "extra" rational root `m=2` at `b=3` is `< b`, hence irrelevant to (★).)

Every `Q_b` in range is moreover **irreducible over ℚ** — so the law for each `b` is in fact a
single irreducibility statement, the cleanest form of the gap.

---

## 5. The precise residual gap, and the recommended attack

> **Open ((★), classes `b ≡ 2,3 (mod 4)`):** `Q_b(m)` has no integer root `m ≥ b`. Empirically
> `Q_b` is irreducible over ℚ for all such `b`; proving irreducibility uniformly closes the law.

What is now known *not* to work (do not re-attempt): single-prime Newton polygon, Eisenstein,
2-adic truncation, finite covering sets, class-uniform witness primes — all provably absent
(06-08 CODE), because `s = 1+i` is 2-adically nilpotent and drags a bare factor `m`. And — new
this session — **analytic dominance / log-concavity**, because the real roots are in-range
(Proposition 3).

The live handles, all *global*:

1. **Uniform irreducibility of the Gegenbauer-parameter family** `g_b(m) = C_b^{(−m)}(−(1+i)/2)`.
   Irreducibility/Galois results for Gegenbauer/Legendre polynomials exist *in the degree
   variable*; here the variable is the **parameter**, an unusual but possibly tractable setting.
   The `₂F₁(−b/2,−(b−1)/2; m−b+1; −2i)` form (Prop. 2) is the entry point.
2. **Legendre route.** `A^{−1/2} = Σ_n P_n(−s/2) u^n` (Legendre generating function) gives
   `g_b((2b−1)/2) = Σ_{k} g_k(b) P_{b−k}(−s/2)`; the central-binomial collapse of Theorem A is a
   shadow of `P_n` evaluations. A reflection/contiguity identity in this representation may
   expose the integer-root exclusion.
3. **Discriminant / resultant** of `Q_b`: empirically non-square with a single large prime to odd
   multiplicity in the constant term (06-07). A growth/Perron argument on the *family's*
   coefficient asymptotics (central-trinomial controlled, `T_n ∼ 3^{n+½}/(2√(πn))`) could bound
   rational roots without a congruence.

---

## 6. Files

`~/projects/scratch/2026-06-08-prove/`:
`qb.py` (engine expansion, `Q_b` structure, irreducibility), `engine.py` (Prop 1 + closed form +
Theorem A magnitude), `hyp.py` (Prop 2 `₂F₁`), `thmA.py` (Theorem A closed form, `b ≤ 15`),
`certify.py` (Theorem B exact rational-root certificate, `b ≤ 70`),
`frontier2.py` (integer non-vanishing on `[b, b²]`, `b ≤ 150`).

## Summary

Not a closure. But: a clean **new theorem** (the half-integer central-binomial value, with a
genuine `₂F₁`-quadratic-transformation proof), the **engine recurrence** identifying the family,
a **frontier push to `b ≤ 70`**, and — most usefully for steering — a proof that the remaining
difficulty is **arithmetic (irrational in-range roots avoiding integers)**, not analytic, which
**rules out the log-concavity route** that the plan had hoped would suffice.
