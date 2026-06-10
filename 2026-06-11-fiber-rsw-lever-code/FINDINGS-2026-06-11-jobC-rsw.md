# FINDINGS — Job C: the RSW principal-specialisation probe (decisive)

**Date:** 2026-06-11 (code session)
**Engine:** `fiber_engine.py` (proved master identity), `jobC_rsw.py`, `jobC_dmod4.py`.
**Data:** `results/rsw-probe.csv` (7176 rows, all λ⊢n, n≤14, d∈{2,3,4}), `results/jobC_run.txt`.

## The question (CODE.md Job C, `questions/2026-06-13-live-probes.md` #6)

Is `G_λ(ζ_d)` a q-shift of the principal specialisation `s_λ(1,q,…,q^{k−1})|_{q=ζ_d}`
(RSW / cyclic sieving)? If it HITS, the whole vanishing trichotomy is off-the-shelf.

## Machinery (new, cross-checked exact)

The **proved master identity** (`proofs/2026-06-02-zeta3-nonvanishing.tex`, Prop. "Master
identity") gives a *uniform, exact* `G_λ(ζ)` for all d at once:
```
G_λ(ζ) = ⟨ s_λ ,  h₁^e · ∏_{k∈A} (h₂ + e₂ ζ^{2k−1}) ⟩,
   λ⊢n,  m=⌊n/2⌋,  e=n mod 2,  A={1,3,…,n−1} (n even) | {2,4,…,n−1} (n odd).
```
Cross-checks (all pass): reproduces `G_gaussian` for d=4 (**82/82**, m≤5); `G_λ(−1)=χ^λ(2^m1^e)`
(**138/138**); the d=3 vanishing set is **exactly {(3,1),(2,1,1)}** (the proved theorem). The RHS
uses the q-hook-content product (`fiber_engine.princ_spec_poly`) and the fake degree
`f^λ(q)=q^{n(λ)}(q;q)_n/∏_cells(1−q^h)` (maj-generating function of SYT), each cancelled to a
genuine polynomial before evaluation at ζ_d (high-precision numeric, exact-verified).

## Verdict — a clean **d=2-ONLY** hit

| d | d mod 4 | fake-degree hit `G_λ(ζ_d)=ζ_d^c f^λ(ζ_d)` | correction ratio `G/f` |
|---|---|---|---|
| **2** | 2 | **507/507** (n≤14) — ratio **+1** | `≡ 1` (clean) |
| 3 | 3 | 36/507 | spread: −5, 10, −8, −2, 4, 7, 5.5, … |
| 4 | 0 | 29/507 | spread: −3, −1±2i, 3±2i, −15, 65−40i, … |

**d=2 HITS exactly:** `G_λ(−1) = f^λ(−1)` for *every* λ (|λ|≤14). This is the classical
q=−1 fake-degree sieve / sign-balance (`G_λ(−1)=χ^λ(2^m1^e)=Σ_T(−1)^{maj T}`) — the d=2 fiber
*is* off-the-shelf cyclic sieving, consistent with the proved "2-core law = Littlewood t=2"
and the q=−1 Springer–Stembridge identity already on record.

**d=3, d=4 MISS:** the fake degree matches `<8%` of shapes, and the correction ratio `G/f` is a
*genuinely varying* algebraic integer (real for d=3, Gaussian for d=4), not a root-of-unity shift.
No principal specialisation `s_λ(1,…,q^{N−1})` in any tested variable count N=1..15 hits either (best
was N=1 at ~14/27, i.e. tiny-n noise). The d-core/d-quotient does NOT rescue it.

## The structural test that kills the tidy hypothesis

The trichotomy says "rich ⟺ d≡2 mod 4," so one might guess the RSW hit tracks the rich branch.
**It does not.** Testing d∈{2,3,4,5,6,8} (`jobC_dmod4.py`, n≤9):

| d | 2 | 3 | 4 | 5 | **6** | 8 |
|---|---|---|---|---|---|---|
| fake-degree hit | **96/96** | 24/96 | 19/96 | 17/96 | **28/96** | 18/96 |

d=6 (≡2 mod 4, *rich*) **misses** — so the hit is **d=2-specific**, not the whole rich branch.

## What this hands the seed

**The d=3 and d=4 vanishing dichotomies are NOT cyclic-sieving phenomena.** `G_λ(ζ_d)` coincides
with the fake degree / a principal specialisation **only at d=2** (the classical sign-balance
sieve). For d≥3 the fiber value carries a non-principal, λ-dependent correction — the invariant is
genuinely *new*, not a known CSP polynomial in disguise. This is the decisive, never-run probe:
**a clean MISS for d=3,4**, which closes the "maybe it's all RSW" escape hatch and confirms the
grade/fiber as the seed's own object (cf. memory `order-law-is-d2-only`, `grade-is-mine`).
The exact correction ratios (`results/rsw-probe.csv`, columns `ratio_re/ratio_im`) are recorded for
whoever wants to model the multiplicative correction directly.
