# FINDINGS — odd-prime / Newton-polygon hunt for `b ≡ 2,3 (mod 4)` (Job 2)

**Date:** 2026-06-08 (code session)
**Scripts:** `job2_oddprime_newton.py` (coeffs, no-root primes, disc),
`job2_newton_recheck.py` (CORRECTED Newton criterion),
`job2_covering_periodicity.py` (covering sets + periodicity), shared `dfour_tworow.py`.

For `b ≡ 2,3 (mod 4)` the 2-adic certificate is unavailable (`q_b ≡ m·(m²+m+1)^k mod 2`
has the root `m≡0`). We hunted for an odd-prime obstruction to a rational (= integer,
since `q_b` monic) root of `q_b = primitive(Q_b)`, which is **irreducible** of degree
`⌊b/2⌋` for these classes (cross-checked).

## A correction to the Newton criterion (important)

The right single-prime test for "no rational root" is: **every edge slope of the `p`-adic
Newton polygon is non-integral** (a rational root is a `p`-adic root of integer valuation =
an integer edge slope). The plan's "all `|slope| < 1`" flag is **wrong** — slope `0` is an
integer `< 1` and corresponds to unit roots, which can be rational. Re-running with the
correct *non-integral-slope* test:

> **VERDICT (Newton, corrected).** For every `b ≡ 2,3 (mod 4)`, `6 ≤ b ≤ 42`, and every
> prime `p < 200`, **no** all-non-integral-slope obstruction exists — for the monic `q_b`,
> for its **reciprocal** `m^d q_b(1/m)`, **or** for the **scaled** `q_b(p·m)`. The
> single-prime Newton-polygon route is **dead**, even after breaking monicity. (The earlier
> "all `|slope|<1`" hits were the spurious slope-0 artifact; they certify nothing.)

This hardens last cycle's verdict and closes the "break-monicity" suggestion negatively.

## The only working local obstruction: "no root mod p" — works per-`b`, never uniformly

Since `q_b` is monic, a rational root is an integer root, which reduces mod every `p`. So
`q_b` root-free mod some `p` ⟹ no rational root. Because `q_b` is **irreducible of degree
`≥ 2`**, Jordan's theorem guarantees its Galois group contains a derangement, so infinitely
many `p` give no root mod `p` — **a witnessing prime exists for every individual `b`**
(smallest one is `≤ 100` for all tested `b`; e.g. `b=6→7`, `b=23→41`, `b=42→43`).

**But there is no uniform finite set.** Testing root-freeness of `q_b` mod each candidate
prime `p ∈ {7,…,61}` over the full ranges `b ≤ 62`:

- **No single prime** is a class-uniform witness (confirmed, both classes).
- **No `b`-periodicity:** for no prime `p` is `{b : q_b root-free mod p}` a union of
  arithmetic progressions in `b` (period `≤ 24`). The root-free `b`'s look pseudo-random.
- **No finite covering set:** apparent small covers — `(47,61)` for `b≡2`, `(37,41,43)` for
  `b≡3` — held only because they were fit to the short list `b ≤ 42`. Extended to `b ≤ 62`
  they **evaporate** (e.g. `b = 62` is root-free at *no* candidate prime `≤ 61`); no set of
  size `≤ 4` covers the class. This is the same small-sample coincidence last cycle warned
  about (a degree-`d` poly is root-free mod `p` with prob `≈ 1/e`), now confirmed at scale.

## Discriminant / Galois angle — clean fact, but no uniform lever

`disc(q_b)` is **non-square for all tested `b`** (`b ≤ 34`, both classes; small-prime part
plus the large cofactor's non-squareness), so `Gal(q_b) ⊄ A_d` — extends the record. But:

> **No fixed odd prime divides `disc(q_b)` to odd multiplicity uniformly across either
> class** (checked among primes `< 10⁵`; the discriminants share many small primes but the
> *parities* of their exponents vary with `b`). So there is no clean uniform ramification
> prime to hang a "no linear factor" argument on.

(`disc` factorisations are only partial — the 250–460-digit discriminants were trial-divided
to `10⁵`, leaving large unfactored cofactors; full factorisation is infeasible and
unnecessary, since a *uniform* witness would have to be small.)

## Honest bottom line for the prove phase

- `b ≡ 2,3 (mod 4)` is genuinely **non-local**: single-prime Newton is dead; "no root mod
  p" settles each `b` individually but admits **no** congruence-uniform prime, **no** finite
  Bonciocat covering set, and **no** `b`-periodicity.
- The clean structural facts available: `q_b` irreducible (deg `⌊b/2⌋`), monic, primitive,
  **non-square discriminant** (Galois `⊄ A_d`).
- **Recommended route:** abandon local-prime methods for these classes and pursue a
  **uniform global argument** — either (i) a direct analytic non-vanishing of `I_b(m)` for
  `m ≥ b` (the imaginary-part inequality the law descends from), or (ii) a uniform
  irreducibility / no-linear-factor proof using the family structure of `Q_b`
  (`Q_b(m) = [u^b]`-coefficient growth) rather than a single congruence.

## Data / reproducibility

- `job2_oddprime_newton.py` — Job 2.1 (coeff factorisations), 2.2 (no-root primes), 2.4
  (disc). Uses bounded trial-division (`< 10⁵`) so it never hangs on huge constants.
  *Its built-in `all-|slope|<1` flag is superseded by the corrected criterion below.*
- `job2_newton_recheck.py` — the **correct** non-integral-slope test (monic/reciprocal/scaled).
- `job2_covering_periodicity.py` — covering-set and `b`-periodicity hunt (factor-free,
  via `primitive(Q_b) mod p`), the source of the "no uniform structure" verdict.
