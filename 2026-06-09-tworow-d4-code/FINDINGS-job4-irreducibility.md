# FINDINGS — Job 4: irreducibility census of Q_b, b ≡ 2,3 (mod 4)

**Date:** 2026-06-09 (code session)
**Scripts:** `job4_irreducibility.py`, `job4b_gegenbauer.py`, `_qbcore.py`, `build_cache.py`
**Data:** `results/job4_irred.csv`, `results/Qb_cache.pkl`

## Headline

**Every Q_b for b ≡ 2,3 (mod 4), 2 ≤ b ≤ 119, is irreducible over ℚ — certified, all 60 values.**

This is the decisive computational result of the session, because:

> For b ≥ 6 in this class, deg Q_b = ⌊b/2⌋ ≥ 3. An irreducible polynomial of
> degree ≥ 2 has **no rational root at all**. Q_b is monic, so "no integer
> root ≥ b" (the law (★)) is *implied* by irreducibility. Hence the two-row
> d=4 law holds for **all b ≡ 2,3 (mod 4), 6 ≤ b ≤ 119** — conditional only on
> the computer-verified irreducibility, not on any analytic/2-adic estimate.

The base cases b = 2, 3 have deg Q_b = 1 (Q_2 = m−2, Q_3 = m−2; root m=2 < b for
b=3, =b for b=2) and are handled directly by the existing low-b argument.

## Method (why mod-p, not factor over ℚ)

SymPy's Zassenhaus factorisation over ℚ is too slow at degree ≈ 60 with the
huge integer coefficients Q_b carries. We instead use the standard fast
certificate:

- **(CERT-IRR)** If Q_b mod p is irreducible (a single factor of full degree)
  for some prime p, then Q_b is irreducible over ℚ.
- Fallback **degree-compatibility**: any ℚ-factor degree must be a subset-sum of
  the mod-p irreducible-factor degrees for *every* good p; intersect over many
  primes. If only {0, deg} survive → irreducible.

Result over primes 3 ≤ p < 3000:
- **9 / 60** values certified by a single irreducible reduction (a Frobenius
  d-cycle, e.g. b=34 irred mod 37, b=46 irred mod 47, b=103 irred mod 109,
  b=114 irred mod 131).
- **51 / 60** certified by multi-prime degree-incompatibility.
- 0 uncertified.

(Initial run to p<120 left 8 large b uncertified — degree-≈59 polynomials need a
larger prime pool for a d-cycle to appear, exactly as Chebotarev predicts;
extending to p<3000 closed all of them.)

## Galois / discriminant side

`disc(Q_b)` is a **perfect square only for b = 2, 3** (where it is the disc of a
linear poly, ≡ 1 by convention). For **all b ≥ 6 the discriminant is a
non-square** ⟹ the Galois group of Q_b is **not contained in A_d** ⟹ it contains
a transposition. Combined with transitivity (irreducibility), this points to a
large Galois group (S_d in the generic case) — the natural ambient for a
uniform irreducibility proof.

## Prop 1 (Gegenbauer-in-parameter) — verified symbolically

`job4b_gegenbauer.py` confirms, as an identity of polynomials in the symbolic
parameter m over ℚ(i):

  **g_b(m) = [u^b](1+su+u²)^m = C_b^{(−m)}( −(1+i)/2 )**,  s = 1+i,

for all b ≤ 14, together with the three-term recurrence
b·g_b = s(m−b+1)g_{b−1} + (2m−b+2)g_{b−2}. Derivation: the Gegenbauer generating
function (1−2xt+t²)^{−α} = Σ C_n^α(x) t^n with x = −s/2, α = −m, t = u gives
(1+su+u²)^m exactly. PROVE may cite this cleanly.

## Verdict feeding PROVE

The law for b ≡ 2,3 is now **equivalent to (and provable via) irreducibility of
Q_b**, which is the stronger and cleaner target. The Galois data (transitive +
contains a transposition) is consistent with an S_d-type group; a uniform proof
should aim at irreducibility, with the Newton-polygon discriminant lever (Job 1)
as the engine.
