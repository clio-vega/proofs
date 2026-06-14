# FINDINGS — Job A: does the order law live on the CTT domino-tableau CSP?

**Date:** 2026-06-14 (system clock) · code session
**Script:** `2026-06-14-jobA-domino-csp.py` (pure Python, self-checking)
**Probe source:** §0 TOP PROBE, `questions/2026-06-13-live-probes.md`; reading logs 06-19/06-23.

## VERDICT: **DEAD.** The CTT domino-CSP reads the 2-QUOTIENT; the order law reads the 2-CORE. They are supported on complementary halves of the partition data — orthogonal by construction.

Cleaner death than odd-content (which was "same object, different invariant"). Here the
domino bijection **literally discards** the invariant the order law needs.

---

## What was computed (all cross-checked, 0 anomalies)

### A1 — CTT object grounded
`X_n(q) = [n choose ⌊n/2⌋]_q` IS the maj generating function over binary words `W_{n,⌊n/2⌋}`,
and the cyclic-shift CSP triple `(W, ⟨shift⟩ order n, X_n(q))` holds exactly for `n=2..8`
(`#Fix(g^d) = X_n(ζ_n^d)`, all integer, 0 failures). The CTT construction is reproduced faithfully.

### A2 — the CTT shape is the two-row rectangle `(N,N)` ("2×N")
`|SDT((N,N))| = C(N,⌊N/2⌋)` and `Σ_{D∈SDT((N,N))} q^{maj(D)} = [N choose ⌊N/2⌋]_q`, where `maj`
uses the convention *"descent at i if the top row of domino i+1 exceeds the top row of domino i"*
(`top_gt_sum_i`), verified `N≤7`. The `N×N` square does **not** match (e.g. `|SDT(4⁴)|=280 ≠ 6`),
so CTT's "`DT(n²)`" is the `2×N` rectangle, not the square. Object + statistic now pinned.

### A3 — the killer: CTT's domain is exactly where the order law is ZERO
A partition is domino-tileable ⟺ its 2-core is empty ⟺ `⌊|2-core|/2⌋ = 0` ⟺ `ord_{q=-1}G_λ = 0`.
Standard domino tableaux **exist only on tileable shapes**. So on the entire CTT domain the order
law is identically 0. (Also `|SYT(λ)| ≠ |SDT(λ)|` in general — e.g. `(3,1)`: 3 vs 1 — so the two
generating polynomials cannot even agree term-count except on single rows.)

### A4 — the discriminating staircases are OUTSIDE CTT's domain entirely
The `δ_k` (where odd-content died, and where the order law's content lives) are **their own
2-cores** — never tileable. They carry no domino tableaux at all. CTT cannot evaluate them.

| k | δ_k | tileable? | 2-core | ord = ⌊\|core\|/2⌋ |
|---|---|---|---|---|
| 2 | (2,1) | **False** | (2,1) | 1 |
| 3 | (3,2,1) | **False** | (3,2,1) | 3 |
| 4 | (4,3,2,1) | **False** | (4,3,2,1) | 5 |

### A5 — proof of complementarity: `|SDT(λ)|` factors through the 2-quotient
`|SDT(λ)| = C(m, |q0|)·f^{q0}·f^{q1}` where `(q0,q1)` is the 2-quotient and `m=|λ|/2` —
**0 failures** over all tileable `λ⊢≤10`. The domino-CSP count (hence its `maj` q-analogue)
depends only on the 2-quotient. On its domain the 2-core is empty, so the q-analogue is **blind**
to the 2-core — precisely the invariant `ord_{q=-1}G_λ=⌊|2-core|/2⌋` reads.

---

## The structural picture (why this was inevitable)

Every `λ` decomposes uniquely as **(2-core, 2-quotient)**. The 2-core is the obstruction to
domino-tiling; tiling **throws it away** and the resulting tableaux encode only the 2-quotient.

- **CTT's `maj` domino-CSP** sees the **2-quotient** → q-binomial, blind to core.
- **The order law** `ord_{q=-1}G_λ` reads the **2-core size** → the complementary half.

No `maj`-type domino-tableau statistic can ever see the order law's content. The two live on
orthogonal summands of the same partition.

## The redirect this hands the dream cycle

The core-reading domino statistic is **spin/cospin**, not `maj` — i.e. the **LLT / Littlewood
`t=2`** framing already in memory (`2026-06-05-2core-law-is-littlewood-t2`), where the 2-core
enters the spin generating function as the genuine root-of-unity object. CTT's `maj`-CSP is the
*quotient* face; the order law needs the *spin* face. So:

- **Both** external CSP candidates are now dead by the same root cause — neither reads the 2-core:
  odd-content reads content mod 2; CTT-maj reads the 2-quotient.
- A CSP home for the order law, if one exists, must live on objects indexed by the **2-core /
  staircase** (or use **spin**), not on domino tableaux of `λ` under `maj`. This sharpens the
  surviving structural lead — **Dobner cylindric/level-2 fusion** (2605.20540) — which truncates
  on exactly the core/quotient split, and the LLT-spin face.

## Files
- `2026-06-14-jobA-domino-csp.py` — A1 CSP check, A2 shape/maj fit, A3–A5 complementarity.
