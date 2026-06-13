# FINDINGS — Job A: the odd-content probe (2026-06-13 code)

**Script:** `2026-06-13-jobA-oddcontent.py` · pure Python (no SageMath in container).

## Question

Armon–Swanson: `f^λ(q,t) = f^λ(q)·∏_{□∈λ}(1 + t·q^{c(□)})`. At `q=−1`, the factor
`(1+t·(−1)^{c})` is `1−t` when `c` is **odd**, `1+t` when even; so the order of the content
product at `(q,t)=(−1,1)` is exactly the **odd-content count**
`OC(λ) := #{□ : c(□) odd}`.

**Test:** is the proved grade order law `ord_{q=−1} G_λ = ⌊|2-core(λ)|/2⌋` a clean function of
`OC(λ)` or `OC(core₂(λ))`? A match would give a 4th, rep-theoretic (Pfannerer-CSP) proof of the
order law, and explain its **d=2-onlyness** (only the 2nd root reads content mod 2).

## Verdict: **MISMATCH — decisively.** The order law is NOT the odd-content count.

### Cross-check first (the machinery is sound)

`ord_{q=−1} G_λ` computed directly from `G_λ(q)=Σ_{T∈SYT(λ)} q^{s(T)}` (parity-twisted descent
statistic `s`) agrees with the closed form `⌊|2-core|/2⌋` on **all 0 failures, n ≤ 9**. The
order law itself is reconfirmed; the mismatch below is real, not a coding artefact.

### The staircase tell (the cleanest possible match — and it fails)

For the 2-core `δ_k=(k,k−1,…,1)`, the contents run `−(k−1)…(k−1)`:

| k | \|δ_k\| | `ord=⌊\|δ_k\|/2⌋` | `OC(δ_k)` | match |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | ✓ |
| 1 | 1 | 0 | 0 | ✓ |
| 2 | 3 | **1** | **2** | ✗ |
| 3 | 6 | **3** | **2** | ✗ |
| 4 | 10 | **5** | **6** | ✗ |
| 5 | 15 | **7** | **6** | ✗ |
| 6 | 21 | **10** | **12** | ✗ |
| 7 | 28 | **14** | **12** | ✗ |
| 8 | 36 | **18** | **20** | ✗ |

Closed forms (worked out and confirmed):
- `OC(δ_k) = ⌊k/2⌋·(⌊k/2⌋+1)`  (= `2·T_{⌊k/2⌋}`, twice a triangular number).
- `ord(δ_k) = ⌊k(k+1)/4⌋`.

These first disagree at **δ₃ (2 vs 3)** and never reconcile. `ord = ⌊|core|/2⌋` is simply half the
core's *size* — it has nothing to do with the even/odd content split. (`min(OC,EC)`, `EC`, etc. all
fail too: e.g. δ₃ has `EC=4`, `OC=2`, neither is 3.)

### Full census (n ≤ 14, 507 partitions, closed form for ord)

- `ord ≠ OC(λ)`  : **504 / 507** mismatches (OC(λ) = OC(core) + #dominoes, grows with n).
- `ord ≠ OC(core₂(λ))` : **120 / 507** mismatches — systematic, governed by the staircase
  formula above. So even restricting to the 2-core, the odd-content count is a *different* clean
  function of the core than `⌊|core|/2⌋`.

## What this kills / what survives

- **KILLS** the hope that the smaj content-product reads off the grade's `q=−1` order — i.e. the
  "super-maj content-product = grade home as a route to the order law" reading of
  `[[thmB-triple-convergence]]`. The factor `∏(1+t q^{c})` at `q=−1` computes `OC`, and `OC ≠
  ⌊|core|/2⌋`. No rep-theoretic proof of the **order law** comes for free this way.
- **Does NOT touch** the multiset identity `{d_j} = {s(T)}` itself (a statement about the whole
  grade multiset, not its `q=−1` order). That leg of the convergence is untested here.
- **Mechanistic note:** the reason is structural. `q=−1` on the *content product* isolates content
  **mod 2** (→ OC). The grade order law isolates the **2-core size** (a domino/abacus invariant).
  These coincide only at the trivial cores δ₀,δ₁. The "d=2-onlyness" of the order law is genuine,
  but it is NOT the content-mod-2 onlyness of the Armon–Swanson product — they are different
  d=2 phenomena that happen to share the prime 2.

**Bottom line:** falsifiable, falsified. The odd-content count is not a back door to the order law.
