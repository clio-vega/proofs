# The order law as a combinatorial minimum: `min_{T in SYT(λ)} s(T) = τ(τ+1)/2`

**Clio — 2026-05-31 (prove session)**

## Theorem (A)

Fix `n >= 1` and a partition `λ ⊢ n` with `ℓ = λ'_1` rows. For a standard Young tableau
`T ∈ SYT(λ)` write `Des(T) = { i ∈ {1,…,n-1} : i+1 lies in a strictly lower row of T than i }`,
and define the **parity-twisted descent statistic**

```
   s(T) = Σ_{i ∈ Des(T)} w_i ,        w_i = 2i-1  if  n-i is odd,   else  w_i = 0 .
```

Let `τ = max(0, 2ℓ - n - 1)`. Then

```
   min_{T ∈ SYT(λ)}  s(T)  =  τ(τ+1)/2 .
```

Equivalently (since `ord_{x=q²} Z_λ = min_j d_j = min_T s(T)`, established earlier) this is the
order law `ord Z_λ = τ(τ+1)/2`, now a statement about tableaux — no representation theory
enters the proof of the lower bound or the `τ>0` construction below.

Call a position `i ∈ {1,…,n-1}` **free** if `n-i` is even (so `i ≡ n mod 2`, `w_i = 0`) and
**paid** if `n-i` is odd (`w_i = 2i-1 > 0`). The paid positions in increasing order
`p_1 < p_2 < ⋯` are: `n` even ⇒ paid `= {1,3,5,…}`, `p_t = 2t-1`, `w_{p_t} = 4t-3`; `n` odd ⇒
paid `= {2,4,…}`, `p_t = 2t`, `w_{p_t} = 4t-1`. In both cases `w_{p_1} < w_{p_2} < ⋯`.

---

## 1. The prefix-descent lower bound (Lemma A)

> **Lemma A.** For every `T ∈ SYT(λ)` and every `m ∈ {1,…,n-1}`,
> `|Des(T) ∩ {1,…,m}| >= r(m) - 1`, where `r(m) = min{ r : λ_1+⋯+λ_r >= m+1 }` is the least
> number of rows that can hold `m+1` boxes.

**Proof.** In any SYT the cells containing `1,…,k` form a Young diagram for each `k` (the cells
with entry `<= k` are closed under moving up or left). Apply with `k = m+1`: entries `1,…,m+1`
occupy a Young subdiagram `μ ⊆ λ`, `|μ| = m+1`. Let `h` be its number of rows. If
`h <= r(m)-1` then `m+1 = |μ| = Σ_{i=1}^h μ_i <= Σ_{i=1}^{r(m)-1} λ_i <= m` (last step by
definition of `r(m)`), contradiction; so `h >= r(m)`.

Insert `1,…,m+1` in order; say `i+1` *opens a row* if its row was empty among `1,…,i`. Exactly
`h` rows are opened, entry `1` opens row 1, so `h-1` openings occur at entries `i+1` with
`1 <= i <= m`. If `i+1` opens row `r`, then (since `{1,…,i+1}` is a Young diagram with newest
nonempty row `r`) rows `1,…,r-1` are nonempty in `{1,…,i}` and rows `>= r` empty there; so
`row(i) <= r-1 < r = row(i+1)`, i.e. `i ∈ Des(T)`, `i <= m`. These give `h-1` distinct descents
in `{1,…,m}`, so `|Des(T) ∩ {1,…,m}| >= h-1 >= r(m)-1`. ∎

Consequences (at `m = n-1`, where `r(n-1) = ℓ`): `|Des(T)| >= ℓ-1`; and `D(m) := r(m)-1`
satisfies `D(n-1) = ℓ-1`. (`D(m) = r(m)-1` is exactly `min_T |Des(T) ∩ [1,m]|`, brute-verified
`n <= 9` in `lowerbound.py`.)

---

## 2. The lower bound `s(T) >= τ(τ+1)/2`

Let `f(m) = #{ i : 1 <= i <= m, i free }`.

> **Lemma B (forced paid descents).** Every `T ∈ SYT(λ)` has at least `⌈τ/2⌉` paid descents.

**Proof.** By Lemma A at `m=n-1`, `|Des(T)| >= ℓ-1`; at most `f(n-1)` of these are free, so the
paid descents number at least `(ℓ-1) - f(n-1)`. Now `n` even ⇒ `f(n-1) = n/2 - 1`, giving
`(ℓ-1)-f(n-1) = (2ℓ-n)/2 = (τ_raw+1)/2` with `τ_raw := 2ℓ-n-1`; `n` odd ⇒ `f(n-1)=(n-1)/2`,
giving `(2ℓ-n-1)/2 = τ_raw/2`. When `τ_raw>0` (`τ=τ_raw`), since `τ_raw ≡ n-1 (mod 2)` this is
`⌈τ/2⌉`. When `τ_raw<=0` (`τ=0`), `⌈τ/2⌉=0` and the claim is vacuous. ∎

> **Lemma C (arithmetic).** `Σ_{t=1}^{⌈τ/2⌉} w_{p_t} = τ(τ+1)/2`.

**Proof.** For `τ>0`: `n` even ⇒ `τ=2k-1`, `Σ_{t=1}^k (4t-3) = k(2k-1) = τ(τ+1)/2`; `n` odd ⇒
`τ=2k`, `Σ_{t=1}^k (4t-1) = k(2k+1) = τ(τ+1)/2`. For `τ=0`, both sides `0`. ∎

> **Proposition (Lower bound).** `s(T) >= τ(τ+1)/2` for every `T ∈ SYT(λ)`.

**Proof.** Let `P = {i ∈ Des(T) : i paid} = {i_1 < ⋯ < i_{|P|}}`. Each `i_t` is a paid position
and `i_1,…,i_t` are `t` distinct paid positions, so `i_t >= p_t`, whence `w_{i_t} >= w_{p_t}`.
With `K = ⌈τ/2⌉ <= |P|` (Lemma B) and weights `>= 0`,
`s(T) = Σ_{t<=|P|} w_{i_t} >= Σ_{t<=|P|} w_{p_t} >= Σ_{t<=K} w_{p_t} = τ(τ+1)/2` (Lemma C). ∎

The lower bound uses Lemma A **only at `m = n-1`**. Independently, `capstone.py` brute-verifies
`min over SYT = LP-bound = τ(τ+1)/2` for every partition of every `n <= 9`.

---

## 3. Achievability `s(T*) = τ(τ+1)/2`

It suffices to exhibit `T* ∈ SYT(λ)` with **no descent at any paid position `> τ`**: then the
paid descents lie in `{p_1,…,p_K}` (paid positions `<= τ`, as `p_K = τ`), so
`s(T*) <= Σ_{t<=K} w_{p_t} = τ(τ+1)/2`, and the lower bound forces equality.

### 3a. Degenerate cases
`ℓ=1`: only SYT has `Des=∅`, `s=0=τ(τ+1)/2`. `λ=(1ⁿ)`: only SYT has `Des={1,…,n-1}`, `τ=n-1`,
every paid position is `<= n-1 = τ`, so `s = Σ_{paid} w = τ(τ+1)/2` (Lemma C). Both match.

### 3b. Explicit construction for `τ > 0` (fully proved)
Here `2 <= τ+2 = 2ℓ-n+1 <= ℓ`. Fill `λ`:
1. **Spine:** `1,…,τ+2` down column 1, rows `1,…,τ+2`.
2. **Alternation:** "up" entries `τ+3, τ+5, …` fill columns `>= 2` in column-major order
   `(1,2),(2,2),…,(λ'_2,2),(1,3),…`; "down" entries `τ+4, τ+6, …` continue column 1 in rows
   `τ+3,…,ℓ`. Counts match (`n-ℓ` up cells, `n-ℓ-1` down cells; tail length `2(n-ℓ)-1`).

*Valid SYT:* column 1 increases down; columns `>= 2` column-major ⇒ increase down and right;
`entry(r,1) < entry(r,2)` (for `r<=τ+2`, left entry `<= τ+2 <` any tail entry; for `r>=τ+3`,
`2r-τ-2 < 2r+τ+1`). *Descents:* the `j`-th up entry has row `r_j <= τ+1+j` (column 2: `r_j=j`;
later columns: `r_j <= τ+1+j`), while the preceding box is in row `τ+1+j`; so each up move is a
non-descent and each down move strictly descends. Hence
`Des(T*) = {1,…,τ} ∪ {τ+1, τ+3, …, n-2}`. Every element `> τ` is `τ+1, τ+3, …`, all **free**
(`τ+1 = 2ℓ-n ≡ n mod 2`, step 2). So no paid descent `> τ` and `s(T*) = τ(τ+1)/2`. ∎

This construction is also computer-checked: `cons3.py` confirms `s = target` for every `τ>0`
shape with `n <= 10`.

### 3c. The remaining (`τ = 0`) shapes — CLOSED (2026-06-02, descent-set Kostka number)

For `τ = 0` (i.e. `2ℓ <= n+1`) the target is `s(T*) = 0`: an SYT with **all descents free**,
i.e. `Des(T) ⊆ S` where `S = { i ∈ [1,n-1] : i ≡ n (mod 2) }` is the set of free positions.
This is now **PROVED uniformly**, not by an explicit construction but by an existence/counting
argument that is cleaner than any construction:

> **Theorem (τ=0 achievability).** If `τ(λ)=0` then `#{ T ∈ SYT(λ) : s(T)=0 } = K_{λ,μ} > 0`,
> where `μ = (2^{n/2})` for `n` even and `μ = (2^{(n-1)/2}, 1)` for `n` odd, and `K_{λ,μ}` is
> the Kostka number (number of SSYT of shape `λ`, content `μ`).

**Proof.** By the classical descent-set identity (Stanley, *EC2*, Cor. 7.19.5 / Thm 7.19.7;
Gessel): for any `S ⊆ [n-1]`, `#{T ∈ SYT(λ) : Des(T) ⊆ S} = K_{λ, α(S)}`, where `α(S)` is the
composition of `n` with partial-sum set `S`. Our `S` (free positions) gives
`α(S) = (2,2,…,2)` for `n` even and `(1,2,2,…,2)` for `n` odd; rearranged to the partition `μ`
above. Now `K_{λ,μ} > 0 ⟺ λ ⊵ μ ⟺ μ' ⊵ λ'`. The conjugate `μ' = (⌈n/2⌉, ⌊n/2⌋)` has only two
parts, so the **only** binding dominance inequality is the first: `μ'_1 = ⌈n/2⌉ ≥ λ'_1 = ℓ`,
i.e. `2ℓ ≤ n+1`, i.e. `τ=0`. Hence `τ=0 ⟹ K_{λ,μ} ≥ 1`, giving an SYT with `s(T)=0`. ∎

So the τ=0 case is not merely "exists"; the **number** of minimizers is exactly the Kostka
number `K_{λ,μ}`. The earlier "open gap" — the greedy that failed 7 wide shapes (smallest
`(3,2,2)`) — was a construction artifact: existence never needed a greedy, only Kostka
positivity, which is pure dominance. Conjugation `s(T') = C(n,2) - s(T)` gives the dual count:
`#{ T' ∈ SYT(λ') : every paid position a descent } = K_{λ,μ}` as well.

**Correction of attribution.** This gap was previously conjectured to be Swanson's
modular-major-index existence theorem (arXiv:1701.04963). On reading, that is *not* the tool:
Swanson controls `maj(T) mod n` — a single linear functional reduced mod `n` — whereas `s(T)=0`
is a constraint on *which positions* (by parity) carry descents, a strictly finer condition on
the descent set. Swanson is one functional too coarse. The correct tool is the elementary
descent-set Kostka positivity above. Verified: `#{T:s(T)=0} = K_{λ,μ}` (0 mismatches, all 194
shapes `n≤11`); `λ⊵μ ⟺ τ=0` (0 failures, all 507 shapes `n≤14`); single binding constraint
`ℓ≤⌈n/2⌉` (0 failures, `n≤12`) — `verify_kostka_bridge.py`. Full writeup:
`2026-06-02-tau0-achievability-kostka.tex`.

---

## 4. Computational verification
`scratch/2026-05-31-minimizer/`: `capstone.py` (min over SYT = LP bound = `τ(τ+1)/2`, `n<=9`),
`lowerbound.py` (`D(m)=r(m)-1`, `n<=9`), `cons3.py` (lower-bound certificate for all SYT, and
the §3b construction; passes every `τ>0` shape, exposes the 7 `τ=0` greedy failures), `n<=10`.
Lemma C arithmetic separately checked to `τ=59`.

---

## 5. Corollary: `τ` is an output of the bound
The lower bound used only `D(n-1)=ℓ-1` and `f(n-1)`; the forced paid count came out as
`(ℓ-1)-f(n-1) = ⌈(2ℓ-n-1)/2⌉`, and the minimum is governed by
`τ = max(0, 2λ'_1 - n - 1)` — re-deriving the column-1 charge formula from prefix congestion
alone, independently of the spectral definition of `τ`.

---

## 6. Status
* **Lower bound `s(T) >= τ(τ+1)/2`: complete and rigorous** (Lemmas A–C + Proposition), via
  Lemma A at the single threshold `m=n-1`; certificate brute-verified `n<=9`.
* **Achievability for `τ > 0`: complete and rigorous** (explicit construction §3b, proved and
  computer-checked `n<=10`).
* **Achievability for `τ = 0`: CLOSED (2026-06-02).** `#{T : s(T)=0} = K_{λ,μ} > 0` for `τ=0`
  via the descent-set Kostka number and the two-part conjugate dominance `ℓ ≤ ⌈n/2⌉` (§3c).
  No construction needed — pure Kostka positivity. **Theorem (A) is now complete and uniform
  over all shapes**: lower bound (Lemmas A–C) + `τ>0` construction (§3b) + `τ=0` Kostka (§3c).

**The order law `ord_{x=q²} Z_λ = τ(τ+1)/2` is now fully proved combinatorially** — no
representation theory, uniform over all partitions, bypassing the merged-arc reach gap that
stalled the operator/deletion-bridge route ~15 sessions.

The companion statement — the full multiset bridge `{d_j}(λ) = { s(T) : T∈SYT(λ) }`, of which
this theorem is the `min`-shadow — remains the prize (Theorem (B) in PROVE.md).
