# FINDINGS — closing `b ≡ 0 (mod 4)` for the two-row d=4 law (Job 1.1 + 1.3)

**Date:** 2026-06-08 (code session)
**Scripts:** `job1_mod4_reduction.py` (Job 1.1), `job1_v2_digits.py` (Job 1.3),
shared `dfour_tworow.py`.
**Companion:** `FINDINGS-fp2-recurrence.md` (Job 1.2, the 𝔽₂ closed form).

Notation as in the companion: `I_b(m) = Im[u^b]((1−u)(1+su+u²)^m)`, `s=1+i`,
`R = ⌊(b−1)/2⌋`. For `b` even, `R = b/2 − 1`.

---

## Job 1.1 — the mod-4 reduction is EXACT, and the reason is `b`-uniform

**Claim verified.** For `b ≡ 0 (mod 4)`, as a congruence of integer-valued polynomials
(all Mahler coefficients agree mod 4 — a rigorous finite certificate):
```
I_b(m) ≡ −m·C(m−1,R) + 2·C(m,2)·C(m−2,R) − 2·C(m,3)·C(m−3,R−1)   (mod 4).
```
Verified **EXACT** for every `b ∈ {8,12,16,…,48}` (degrees of the difference up to 47).

**Mechanism (why `j ≥ 5` and `j ≡ 0 mod 4` terms vanish mod 4).**
`I_b` is a signed sum over `j` of `C(m;j,c)·Im((1+i)^j)` with `c` fixed by `j+2c ∈ {b,b−1}`.
Now `Im((1+i)^j) = 2^{j/2} sin(jπ/4)`, whose value mod 4 is:

| `j`               | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|-------------------|---|---|---|---|---|---|---|---|---|---|----|----|----|
| `Im((1+i)^j)`     | 0 | 1 | 2 | 2 | 0 |−4 |−8 |−8 | 0 |16 | 32 | 32 | 0  |
| `mod 4`           | 0 | 1 | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0  | 0  | 0  |

Only `j ∈ {1,2,3}` survive mod 4, **independently of `b`**. They give exactly the three
displayed terms (`j=2` from the `1`-part, `j=1,3` from the `−u`-part). So the reduction
holds for **all** `b ≡ 0 (mod 4)`, not merely `b ≤ 48` — the finite check confirms a uniform
mechanism. `v₂(Im((1+i)^j)) ≥ 2` for every `j ≥ 4` with non-zero imaginary part
(`v₂ = ⌊j/2⌋`), which is the clean statement the prove phase can cite.

---

## Job 1.3 — the valuation proposition holds, and `b ≡ 0 (mod 4)` essentially closes

**Proposition verified.** For `b ∈ {8,12,16,20,24}` and **every** `m ∈ [b, b+400]`:
```
v₂( I_b(m) )  =  Σ_{r=0}^{R} v₂(m−r)  −  v₂(R!)     ( = v₂((m)_{R+1}) − v₂(R!) ),
```
i.e. `Q_b(m) = I_b/(m)_{R+1}` has constant 2-adic valuation `−v₂(R!)` — a 2-adic **unit**
times `1/R!`. `δ := lhs − rhs ≡ 0` throughout, with **no zeros of `I_b`** in any window.

**Which `j`-terms are valuation-minimal (the tie `τ₁`) — UNIFORM in `b`.**
Decomposing `I_b(m) = Σ_j D_j(m)` (each `D_j` an exact integer) and finding the
minimal-`v₂` terms, across **all** `b ∈ {8,…,48}` and `m ∈ [b,b+400]`:

| parity of `m` | valuation-minimal `j`-set | count of minimisers |
|---------------|---------------------------|---------------------|
| `m` even      | `{1}`                     | **1** (unique)      |
| `m` odd       | `{1, 2, 3}`               | **3**               |

This is the exact, `b`-uniform input the prove phase needs. It yields the **proof skeleton**:

- **`m` even:** `j=1` (the term `−m·C(m−1,R)`) is the *unique* valuation-minimiser, so
  `I_b(m)/2^{v_min}` is odd ⟹ `I_b(m) ≠ 0`.
- **`m` odd:** exactly `j ∈ {1,2,3}` tie at `v_min` and all other terms are strictly higher,
  so `I_b(m)/2^{v_min} ≡ (D_1+D_2+D_3)/2^{v_min} ≡ (odd)+(odd)+(odd) ≡ odd (mod 2)` —
  **three odd numbers sum to an odd number** ⟹ a 2-adic unit ⟹ `I_b(m) ≠ 0`.

So for `b ≡ 0 (mod 4)` the entire order law reduces to the single **`b`-uniform**
combinatorial statement *"the valuation-minimal index set is `{1}` (`m` even) or `{1,2,3}`
(`m` odd)"*, verified here for `b ≤ 48` over windows of length 400. The remaining task for
the prove phase is to establish that argmin structure uniformly in `(b,m)` — a finite,
explicit 2-adic computation on three binomial valuations `v₂(C(m−1,R)), v₂(C(m−2,R)),
v₂(C(m−3,R−1))` (the "mod-4 lemma" of the memory record), now pinned down precisely.

---

## Bottom line

- Job 1.1: the mod-4 reduction is **exact and `b`-uniform** (mechanism: only `j∈{1,2,3}`
  survive mod 4).
- Job 1.3: the valuation proposition holds with **no exceptions** over long windows, and the
  argmin tie is **always `{1}` (`m` even) / `{1,2,3}` (`m` odd)** — the "3 odds = odd" unit
  argument closes `b ≡ 0 (mod 4)` modulo a finite uniform binomial-valuation lemma.
- Combined with the already-proved `b ≡ 1 (mod 4)`, this leaves only the genuinely hard
  `b ≡ 2,3 (mod 4)` (Job 2), where 2-adics are provably insufficient.
