# FINDINGS — Jobs A & B: Route-B pruned; the δ_j ladder, anchor, and SYT model

**Date:** 2026-06-11 (code session)
**Engine:** `jobB_ladder.py`, `jobA_confirm.py` (reuse `job1_tie_census.py`, `core_quotient.py`).
**Data:** `results/e2-rung-ladder.csv` (43 255 rows, m≤12).

> **Framing note (honest):** CODE.md asked these jobs to *arm* a proof of the step law M★. This
> cycle's PROVE session proved M★ is **tautological** on J\*-pairs (`val` constant on J\* forces
> `Δv₂(M)=−2^{a−1}−Δv₂(bin)`; it carries no arithmetic content). So the data below is recorded as
> *genuine structural data on `M_j`* — not as a lever for M★. The pure-M 8-step check is included
> only as a tautology cross-check (38/38, flagged as 0=0).

## Job A — Route B (Ayyer–Kumari hook-Schur factorisation): **PRUNED** (confirmed)

Browsing is off, so the AK statement cannot be re-transcribed; PROVE already pruned it
(`steplaw doc §4`). CODE confirms the two structural facts that kill the import:

1. **(2,2) is itself a 4-CORE** — `4-quotient((2,2))=((),(),(),())`. AK factorise/predict
   *non-vanishing* on empty-4-core shapes; (2,2) is empty-4-core yet `G_{(2,2)}(i)=0`. So
   "import via empty 4-core" predicts the wrong answer at *the only vanisher*.
2. **No product structure to carry the −4.** Tool C gives `M_j=Σ_μ K^{(j)}_{λμ} f^μ`, a **positive
   sum** of SYT-counts, not a product. On the 38 pure-M J\*-pairs, `M_{j+8}` is a sum of ≥2 distinct
   `f^μ` (35/38). There is no single linear factor whose `v₂` could account for the jump; the −4 is
   a Newton-locus (J\*-membership) fact, not a factorisation fact. **Route B stays pruned.**

## Job B — the e₂-rung ladder `δ_j = v₂(M_{j+1})−v₂(M_j)`

`results/e2-rung-ladder.csv`: per-rung `δ_j`, `Δv₂(bin)`, J\*-membership, all λ⊢2m, m≤12.
`δ_j` is broad and roughly symmetric about 0 (range −16..+16; mode 0 at 11 332, tails to ±16).
**No per-rung law** — and the SYT model (below) explains why: `M_j` is a positive sum, so `v₂(M_j)`
is `v₂` of a sum, not a min of valuations. This re-confirms Job A (06-10): *no global `v₂(M_j)`
closed form exists*, and the per-rung `δ_j` inherits that.

### The anchor `v₂(M₀)=v₂(f^λ)` — confirmed, provable base case

- `M₀ = f^λ` **exactly: 4114/4114** (m≤12).
- **2-quotient hook-length factorisation** of `v₂(f^λ)` holds on every empty-2-core shape
  (**3131/3131**): `v₂(f^λ)=v₂(multinomial(2m;2|q₀|,2|q₁|))+v₂(f^{q₀})+v₂(f^{q₁})`. This is the
  classical James–Kerber/Macdonald dimension factorisation — a clean, provable base for any
  induction up the `j`-ladder.

### Tool C — the positive SYT-count model (confirmed at scale)

`M_j = Σ_μ (#vertical-2-strip chains λ→μ of length j) · f^μ`, via `e₂^⊥ s_ν = Σ_{ν/μ vert-2-strip}
s_μ`. **Matches `M_j`: 2044/2044 (m≤7).** Every `M_j` is a nonnegative integer combination of
SYT-counts with no cancellation — the structural reason `v₂(M_j)` has no min-of-valuations formula.

## What this hands PROVE

The `M_j` ladder is a positive-sum object with a clean provable **anchor** (`M₀=f^λ`, 2-quotient
factorisation) but **no per-rung valuation law** (δ_j is structureless; M★ is tautological). Route B
is pruned. Any induction up the ladder must therefore ride the **sharp Newton polygon** structure
(H1/H2 in the lever findings), not a rung-by-rung `v₂` identity.
