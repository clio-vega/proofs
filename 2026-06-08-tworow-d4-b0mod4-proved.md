# Two-row d=4 law: a 2-adic proof for all b ≡ 0 (mod 4)

**Date:** 2026-06-08 (prove session)
**Target (from PROVE.md):** close `b ≡ 0 (mod 4)` — prove the identity (V)
`v₂(I_b(m)) = v₂((m)_{R+1}) − v₂(R!)`, which gives `Q_b(m) ≠ 0` and hence the two-row
case of the d=4 fiber-vanishing law `G_{(2m−b,b)}(i) = 0 ⟺ (2,2)`.

**Result.** **Proved.** Combined with last cycle's `b ≡ 1 (mod 4)`, the two-row d=4 law now
holds unconditionally for the **full residue classes `b ≡ 0,1 (mod 4)`** — half of all `b`.
The clean LaTeX proof is `2026-06-08-tworow-d4-b0mod4-proved.tex` (compiles, 4 pp).

The mechanism is new and rather pretty: in the `b ≡ 1` case the `j=1` term *strictly*
dominates every other term in `v₂` (the minimal valuation is attained once). For `b ≡ 0`
that domination **fails** — for odd `m`, the terms `j = 2` and `j = 3` *tie* `τ_1`. The law
is rescued not by domination but by **parity counting**: the number of valuation-minimal
terms is always odd (1 for even `m`, 3 for odd `m`), so their sum keeps an odd leading 2-adic
digit. *Three odd numbers sum to an odd number.*

---

## 0. Setup (proved earlier; not re-derived)

`λ = (2m−b, b)`, `s = 1+i`, `P(u) = (1+su+u²)^m`, `I_b(m) := Im[u^b]((1−u)P) ∈ ℤ`.

- **(E)** `I_b(m) = Σ_{j=1}^b τ_j(m)`, `τ_j = C(m,j)·Im((1+i)^j)·(T(b−j)−T(b−1−j))`,
  `T(l) = C(m−j, l/2)` for `l` even else 0. Exactly one of `T(b−j), T(b−1−j)` is nonzero,
  `= ±C(m−j, β_j)`, `β_j = ⌊(b−j)/2⌋`.
- `Im((1+i)^j) = 2^{j/2} sin(jπ/4) ∈ ℤ`; `v₂ = ⌊j/2⌋` for `j ≢ 0 (mod 4)`, `= 0` (term
  vanishes) for `j ≡ 0 (mod 4)`; **odd iff `j = 1`**.
- Forced roots `{0,…,R}`, `R = ⌊(b−1)/2⌋`; `I_b = (m)_{R+1}·Q_b`; the law is
  **(★) `Q_b(m) ≠ 0` for integer `m ≥ b`**, equivalently `I_b(m) ≠ 0` for `m ≥ b`.

For `b ≡ 0 (mod 4)`, `b = 4t`: **`R = (b−2)/2 = 2t−1` is ODD**, `R+1 = b/2` even. The
oddness of `R` drives the whole proof.

## 1. Main theorem

> **Theorem.** For `b ≡ 0 (mod 4)`, `R = (b−2)/2`, and every integer `m`:
> **(V)** `v₂(I_b(m)) = v₂((m)_{R+1}) − v₂(R!)`.
> Hence `v₂(Q_b(m)) = −v₂(R!)` is constant, so `Q_b(m) ≠ 0`, `(★)` holds, and the two-row
> d=4 law holds for all `b ≡ 0 (mod 4)`.

Combined with `b ≡ 1 (mod 4)` (last cycle, same identity (V)): **law holds for all
`b ≡ 0,1 (mod 4)`.**

## 2. The j=1 term

`b` even ⟹ surviving extraction is `T(b−2) = C(m−1,R)`, so
`τ_1 = −m C(m−1,R) = −(R+1) C(m,R+1) = −(m)_{R+1}/R!`, and
`a := v₂(τ_1) = v₂((m)_{R+1}) − v₂(R!)` = RHS of (V). The theorem is `v₂(I_b) = a`.

## 3. Exact term decomposition

For `τ_j ≠ 0` (`j ≢ 0 mod 4`), `h := ⌊j/2⌋`, `s_j := j + β_j = R+1+h`. Vandermonde
`C(m,j)C(m−j,β_j) = C(m,s_j)C(s_j,j)` + telescoping `C(m,s_j)/C(m,R+1) = ∏_{t=1}^h
(m−R−t)/(R+1+t)`, minus `a`, gives the **exact**

> `v₂(τ_j) − a = D₀(j) + S_j(m)`, `S_j(m) = Σ_{t=1}^h v₂(m−R−t)`,

with `m`-independent offset `D₀(j) = h + v₂(R!) − v₂(j!) − v₂(β_j!)`. Evaluating
(`v₂((2h)!) = v₂((2h+1)!) = h + v₂(h!)`; `β_j = R−h` for `j` odd, `R+1−h` for `j` even):

> **Lemma 1.** `D₀(j) = v₂(C(R,h))` for `j` odd; `D₀(j) = v₂(C(R,h−1)) − v₂(h)` for `j`
> even. *(j even uses `h·C(R,h) = (R+1−h)·C(R,h−1)`.)*

## 4. (P1) No term dips below τ_1

> **Lemma 2 (P1).** `v₂(τ_j) ≥ a` for all `j`, all `m`: i.e. `D₀(j) + S_j ≥ 0`.

- `j` odd: `D₀ = v₂(C(R,h)) ≥ 0`, `S_j ≥ 0`. ✓
- `j` even: `m−R−1, …, m−R−h` are `h` consecutive integers, so their product is divisible
  by `h!` ⟹ `S_j ≥ v₂(h!) ≥ v₂(h)`. And `D₀ = v₂(C(R,h−1)) − v₂(h) ≥ −v₂(h)`. Sum ≥ 0. ✓

## 5. Which terms tie

> **Lemma 3.** Among `j ≥ 2` (`τ_j ≠ 0`), `v₂(τ_j) = a` is possible only for `j ∈ {2,3}`,
> and (since `R` odd) `v₂(τ_2) = v₂(τ_3) = a ⟺ m` odd; for `m` even both `> a`; and every
> `j ≥ 4` has `v₂(τ_j) > a` for all `m`.

A tie needs `D₀(j) + S_j = 0`. Lower-bound `S_j`:
- `j` even, `h ≥ 3`: `S_j ≥ v₂(h!) = v₂(h) + v₂((h−1)!) ≥ v₂(h)+1` (since `(h−1)!` even),
  while `D₀ ≥ −v₂(h)` ⟹ sum ≥ 1: **strict**.
- `j` even, `h = 2` ⟹ `j = 4 ≡ 0 (mod 4)` ⟹ `τ_4 = 0`: excluded.
- `j` odd, `h ≥ 2` (`j ≥ 5`): two consecutive integers among `m−R−t` ⟹ `S_j ≥ 1`, `D₀ ≥ 0`
  ⟹ sum ≥ 1: **strict**.
- `j = 2` (`h=1`): `D₀ = v₂(C(R,0)) = 0`, `S_2 = v₂(m−R−1)`.
- `j = 3` (`h=1`): `D₀ = v₂(C(R,1)) = v₂(R) = 0` (`R` odd), `S_3 = v₂(m−R−1)`.

So `v₂(τ_2)−a = v₂(τ_3)−a = v₂(m−R−1)`. `R` odd ⟹ `m−R−1` odd ⟺ `m` odd. ∎

## 6. Proof of theorem

By Lemma 2, `N(m) := I_b(m)/2^a = Σ_j τ_j/2^a ∈ ℤ`, and mod 2 only tying terms survive;
each `τ_j/2^a` for a tying `j` is **odd**. By Lemma 3:
- `m` even: only `j = 1` ties ⟹ `N(m) ≡ 1 (mod 2)`.
- `m` odd: ties are exactly `{1,2,3}` ⟹ `N(m) ≡ 1+1+1 ≡ 1 (mod 2)` (**three odds**).

Either way `N(m)` odd ⟹ `v₂(I_b(m)) = a` = (V). For `m ≥ b > R`, RHS finite ⟹ `I_b(m) ≠ 0`
⟹ `(★)`. ∎

## 7. Verification

`~/projects/scratch/2026-06-08-prove/verify.py`: the decomposition, Lemma 1 closed forms,
(P1), the tie classification, and (V) checked exactly for all `b ≡ 0 (mod 4)`, `4 ≤ b ≤ 64`,
over `m ∈ [b,b+20]` plus 40 random `m ≤ 10⁶` per `b` — **24,480 term-level checks, 0
problems**. `explore.py` prints per-term valuation tables (`b = 4,8,12,16`).

## 8. Status

- **Proved (new):** two-row d=4 law for all `b ≡ 0 (mod 4)`; with `b ≡ 1` ⟹ all
  `b ≡ 0,1 (mod 4)`.
- **Open:** `b ≡ 2,3 (mod 4)` — 2-adic method structurally dead (`q_b` has a genuine root
  mod 2); needs an odd prime (Filaseta Newton-polygon / uniform irreducibility). Unchanged.
