# FINDINGS — Job 4: Hidaka–Itoh (2403.10817) orientation probe

**Date:** 2026-06-09 (code session)  **Script:** `job4_hidaka_itoh.py`

**Two clarifications first (both verified computationally).**
(1) **My fiber `G_λ` is NOT the fake-degree / q-hook polynomial.** Comparing
`G_λ(i)=⟨s_λ,ψ^m⟩` (`ψ=h₁²+i(1+i)e₂`) against `f^λ(q)=Σ_T q^{maj}` at `q=i` gives only
**9/294 matches** (m≤7). My `G_λ` is a genuinely d=4-specific construction, not a
principal specialization. (2) The fake degree is not Hidaka–Itoh's object either:
`f^λ(ζ₂)` already reaches magnitude 8 (n=8) / 20 (n=10), far outside `{−1,0,1}`.

**The probe is therefore about the two arithmetic GATES, not a literal value set.**

- **Hidaka–Itoh gate:** `n` has **≤ 2 distinct odd prime factors**. For small `n` this is
  almost always true — the first failure is `n = 3·5·7 = 105`; `n=210=2·3·5·7` also fails.
- **My richness gate:** `d ≡ 2 (mod 4)` (the `−1/p₂` factor). A 1-in-4 density condition,
  *independent of how many odd primes divide d*.

**Verdict — the gates are distinct arithmetic, neither sharpens the other.**
They coincide only on the smallest rich cases `d ∈ {2, 2p, 2pq}` (which are both rich and
≤2-odd-prime). But my rich family contains `d = 2·p₁p₂p₃` (e.g. **`d=210`: rich
(≡2 mod4) yet FAILS Hidaka–Itoh** with 3 odd primes), and Hidaka–Itoh's family contains
all `d = 4, 8, 9, 2²·p …` (powers of 2 and prime powers) that are **not** rich. So:

> Hidaka–Itoh controls **value boundedness** of a principal specialization via the *count*
> of odd prime factors; my trichotomy controls **vanishing richness** of `G_λ` via the
> mod-4 class. Different phenomena. Hidaka–Itoh does not bear on the d=4 even-|J\*| crux.

(One un-pruned browse lead, now closed as orthogonal — 30-min orientation only.)
