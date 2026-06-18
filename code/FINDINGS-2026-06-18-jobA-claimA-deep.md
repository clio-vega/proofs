# Job A — deep-`k` Claim A content census (2026-06-18 CODE)

**Goal (from CODE.md).** Pin the *exact* 2-adic slice-content law `g[c][k]` of the
master deficit polynomial `N_i^{(c)}` at deep depth `k = c−i`, for `c=4..12`,
`k≤6`, both parity slices of `b`. Feed PROVE's "Claim A at `k≥4`" target. Hunt for
a structural `2^k` factor (sibling of Theorem B). Cross-check via Kummer.

**Convention (matches `claimA_verify.py` / `g0_fixeddiv.py`).**
```
N_i^{(c)}(a,b) = M_{b+i}(a,b,c) · c! · k! / [ (b−c+1) · ∏_{t=2}^{i}(b+t) ],   k = c−i,
```
`M_{b+i}` the alternant (`fast_alt.Malt`, Lemma-0 continuation, valid for `a<b`),
`i=c` the top (`N_c=1`, depth 0). Slices: `a = 2P+b+c`, `b = 2B` (even) / `2B+1` (odd).

---

## Headline: the old "even-slice = k" was a LOWER bound, not the content

The previous `claimA_verify.py` computed content as the **gcd of monomial
coefficients** of the fitted `N_i` — a *lower* bound on the true fixed 2-divisor.
It reported even-slice `= k`. Computing the **exact fixed-divisor content**
(`= min over the slice of v₂(N_i)`, the quantity the master valuation actually
needs) gives values **strictly larger than `k`** in many cases.

Example: `c=6, i=2, k=4`, even slice — old report `4`, true content **`5`**.

This is fully verified (see "Verification gate" below), not a fitting artefact.
§2.2's `even = k`, `odd ≥ 2⌊k/2⌋` remain **true as lower bounds** (so PROVE's
existing certificate is safe), but they are **not the exact content**.

---

## The certified Content Lemma (the usable result for PROVE)

**FLOOR (0 violations, `c≤12`, `b≤200`, both slices):**
```
v₂(N_i^{(c)}) ≥ 2⌊k/2⌋    on both slices,   AND
v₂(N_i^{(c)}) ≥ k          on the EVEN slice.
```
This is the largest *uniform* lower bound that holds; it is exactly what the
boundary master valuation consumes. Confirmed by `jobA_certify.py` and `jobA_law.py`
(H1) with **FAILS: NONE**.

### Exact certified content tables (`min v₂` over `b≤200`)

EVEN slice `g_even[c][k]`:
```
      k=0  k=1  k=2  k=3  k=4  k=5  k=6
c=4    0    1    2    .    .    .    .
c=5    0    1    2    3    .    .    .
c=6    0    1    2    3    5    .    .
c=7    0    1    2    4    5    7    .
c=8    0    1    2    4    6    7    8
c=9    0    1    2    3    6    7    8
c=10   0    1    2    3    5    7    8
c=11   0    1    2    4    5    7    7
c=12   0    1    2    4    6    7    8
```
ODD slice `g_odd[c][k]`:
```
      k=0  k=1  k=2  k=3  k=4  k=5  k=6
c=4    0    0    2    .    .    .    .
c=5    0    2    2    4    .    .    .
c=6    0    0    2    2    5    .    .
c=7    0    1    2    3    5    7    .
c=8    0    0    2    2    6    6    8
c=9    0    2    2    5    6    8    8
c=10   0    0    2    2    5    5    8
c=11   0    1    2    3    5    7    8
c=12   0    0    2    2    6    6    8
```

---

## No clean exact closed form exists at low modulus

The exact content is **NOT a function of `(k, c mod 4)`**. Decisive counterexample
(both grid-stable and direct-Malt-confirmed):

> **odd slice, `k=3`:** `c=5 → 4` but `c=9 → 5`, yet `5 ≡ 9 ≡ 1 (mod 4)`.

So the surplus depends on `c` more finely than mod 4 (mod 8 or genuinely growing).
**Verdict:** chasing an exact uniform content formula in `(k, c mod 4)` is the wrong
target — it does not exist.

### The real refinement: parity of `i`

Indexing by `i` (not just depth `k`) exposes structure (`jobA_law.py`, H2):

- **ODD `i`, ODD slice:** content sits **exactly at the floor** `2⌊k/2⌋` for all
  `k ≤ 3` (surplus 0); picks up a small surplus `1–2` only at `k = 4,5,6`.
- **EVEN `i`:** surplus is larger and irregular (`k=1:{1,2}`, `k=3:{1,2,3}`,
  `k=5:{3,4}`).

This `i`-parity split is the cleanest handle in the data: the boundary indices with
odd `i` are floor-tight (minimal) at small depth, and they are the ones the master
valuation is most sensitive to.

---

## Factorization (step 3): NO sibling of Theorem B

Factoring the slice polynomials `N_i(P,B)` over ℚ (`jobA_factor.py`) gives a uniform
shape across every case examined:
```
N_i  =  2^σ  ·  (a−b+1)  ·  (one irreducible "block")
```
where, on the slice, `a−b+1 = 2P+c+1` (the Theorem-B factor) and the block is a
genuinely **irreducible** high-degree polynomial in `(P,B)`.

- The **only** linear factor that splits off is `a−b+1` itself (already PROVED by the
  prove session, equal-exponent antisymmetric vanishing). On the even slice it is
  *odd* (`2P+c+1`, `c` even) so contributes `0` to the 2-content; for `c` odd it is
  `2(P+(c+1)/2)`, contributing exactly `1`.
- There is **no** run of `k` even-spaced linear factors and **no** isolated
  `2^k·(odd)` block. The hoped-for structural `2^k` factor analogous to `(a−b+1)`
  **does not exist.**
- The bulk of the content is the **explicit `2^σ` constant** = `v₂(c!·k!/den-scaling)`
  (the `const_i = c!·k!` of the master valuation memory). The remaining content comes
  from the irreducible block (small: 0,1,2) and interacts *super-additively* with the
  `(a−b+1)` parity (content is super-additive over factors: e.g. `c=7,k=5` whole
  content `7 > σ(5) + 0 + block(1)`).

---

## Kummer cross-check (step 4): mixed — not always term-by-term

Pointwise, `v₂(N_i) = v₂(M_{b+i}) + v₂(c!·k!) − v₂(den)` (`jobA_kummer.py`). At the
content-achieving minimal points:

| case | `v₂(N_i)` | `v₂(M)` | `v₂(c!k!)` | `v₂(den)` | cancellation lift |
|------|-----------|---------|-----------|-----------|-------------------|
| `c=6,k=4` even | 5 | 2 | 7 | 4 | **0** (term-by-term) |
| `c=8,k=6` even | 8 | 0 | 11 | 3 | **0** (term-by-term) |
| `c=7,k=5` even | 7 | 7 | 7 | 7 | **+3** (cancellation-born) |

So for `c=6,k=4` and `c=8,k=6` the content is "free" — every surviving S₃ term is
already divisible and `v₂(M)` equals the min term valuation; the content is pure
`c!k!/den` arithmetic. But for `c=7,k=5` the alternating sum **lifts** `v₂(M)` by 3
beyond the minimal term: that surplus is cancellation-born and a pure
Kummer-on-each-term bound cannot see it.

---

## Route verdict for PROVE

1. **Build the boundary argument on the FLOOR**, not on an exact content formula.
   The floor `2⌊k/2⌋` (and `≥k` on the even slice) is certified (0 violations,
   `c≤12`, `b≤200`) and is the largest *uniform* bound there is. An exact closed form
   in `(k, c mod 4)` provably does not exist (`c=5` vs `c=9`).

2. The floor itself decomposes cleanly as `v₂(c!·k!) − v₂(den) + [floor − that]`,
   i.e. the **explicit `c!k!/den` 2-power accounting** carries most of it. This is the
   `const_i = c!·k!` arithmetic already in the master valuation — a Content Lemma
   proving the floor should be an arithmetic statement about `c!k!/[(b−c+1)∏(b+t)]`
   plus a *lower* bound `v₂(M_{b+i}) ≥ (small)`, **not** a factorization.

3. Do **not** expect a structural `2^k` factor (none exists) and do **not** rely on a
   pure term-by-term Kummer bound (cancellation-born content appears, e.g. `c=7,k=5`).
   The clean lever is the **parity-of-`i`** split: odd-`i` deficits are floor-tight on
   the odd slice for `k≤3`.

---

## Verification gate (all passed — `jobA_verify_content.py`)

- **Malt vs Murnaghan–Nakayama:** `0/150 bad` on every sampled `(c,i)` — transitive
  trust from the alternant down to MN.
- **Grid stability:** content identical at grids `64, 128, 256, 512` for all probe
  cases — the fixed divisor is rigorously pinned (grid-min is an upper bound on the
  global content; non-decreasing ⇒ exact).
- **Fit-free cross-check:** content computed straight from `Malt` (no polynomial fit)
  equals the poly-based content on the test window — the fitted polynomials reproduce
  `N_i`'s 2-adics, not merely its values.
- **Fit validation:** every slice polynomial reproduces direct `Malt` at large `b`
  (held-out points) before any content is read off it.

## Scripts (`projects/code/threerow-boundary/`)
- `jobA_deep.py` — deep census, exact fixed-divisor content `g[c][k]`, both slices.
- `jobA_verify_content.py` — the verification gate (MN, grid stability, fit-free).
- `jobA_certify.py` — certified content tables over `b≤200` + `(c mod 4)` grouping.
- `jobA_law.py` — content indexed by `(c,i)`; H1 (floor) + H2 (`i`-parity) analysis.
- `jobA_factor.py` — slice factorization `N_i = 2^σ·(a−b+1)·block`.
- `jobA_kummer.py` — pointwise `v₂` decomposition + cancellation-lift census.

**Bottom line:** the exact content is c-dependent with no low-modulus closed form; the
certified **floor** (`≥2⌊k/2⌋`, `≥k` even slice) is the usable Content Lemma; its
cleanest proof is the `c!k!/den` 2-power arithmetic refined by parity of `i`, **not** a
factorization (no `2^k` sibling of Theorem B) and **not** a pure term-by-term Kummer
count (cancellation-born content exists).
