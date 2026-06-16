# The `c = 3` boundary lemma — CLOSED. Three-row `(a,b,3)` is now a COMPLETE theorem

**Date:** 2026-06-16 (prove session)
**Status:** **Complete.** The three boundary residuals (A), (B), (C) flagged in
`2026-06-15-boundary-lemma-threerow.md` §4.2 are all closed, by a hand proof verified numerically.
Together with the proved interior (`2026-06-14-threerow-c3-Jstar-even.md`,
`2026-06-15-compensation-lemma-B-proved.md`) this upgrades the three-row `d=4` family **`λ=(a,b,3)`**
to a **complete theorem** (interior + boundary), joining the already-complete `c=1` and `c=2`
families. The keystone is a **two-factor factor-in-product principle** (two-factor Lemma F), the
general-`c` `a`-odd-top engine.

Companion to `2026-06-15-boundary-lemma-threerow.md` (which closed `c=1,2` and the `c=3` `a`-even top,
and reduced the rest to (A),(B),(C)).

---

## 0. Setup and what remained

Notation as in the companion note. For `λ=(a,b,3) ⊢ 2m` with `a+b = 2m−3` (so `a,b` have **opposite
parity**), set

> `G_j = ⟨s_λ, e₂^j h₁^{2(m−j)}⟩ = C(m,j) M_j`,  `val(j) = j + 2 v₂(G_j)`,  `V = min_j val(j)`,
> `J* = {j : val(j) = V}`,  `Δ(j) := val(j) − val(0)`,  `P := m−b−3 ≥ 0`.

The interior theorem gives `V = val(0) − θ` with offset/`θ` set by the parity of `a`:

| parity | `j₀ = min J*` | `θ` | boundary requirement `val(b+i) > V` ⟺ |
|---|---|---|---|
| `a` even (`b` odd) | `0` | `0` | `Δ(b+i) > 0` for `i = 1,2,3` |
| `a` odd (`b` even)  | `3` | `3` | `Δ(b+i) > −3` for `i = 1,2,3` |

`M_j` is supported on `0 ≤ j ≤ b+3`, so the boundary indices are `j ∈ {b+1, b+2, b+3}`. We work in
the **box-interior range** where the interior 2-adic box lies in `[0,b]`: `b ≥ 6` (`a` even, so
`b ≥ 7` odd) and `b ≥ 9` (`a` odd, so `b ≥ 10` even). The finitely many smaller `b` are the
per-family boundary-tie shapes already handled directly in the `c=3` interior note.

**The standing relations** (from `a = 2P+b+3`):

> `a − 1 = 2P + b + 2`,  `a − b + 1 = 2(P+2)`,  `a + 2 = 2P + (b+5)`.

In particular `v₂(a−b+1) = 1 + v₂(P+2)`, and **when `b = 2β` is even** (`a` odd),
`a − 1 = 2(P+β+1)`, so `v₂(a−1) = 1 + v₂(P+β+1)`.

### 0.1 The broken assumption

The companion note conjectured it would "suffice" to prove the standalone valuation bounds
`v₂(N₂) ≥ v₂(P+1) − 1` (claim A) and `v₂(N₁) ≥ v₂(P+2) − 1` (claim B), where

> `N₂ = a(b+3) − (b²+2b+3)`,  `N₁ = ab² + 5ab + 6a − b³ − b² − 4b − 6`

are the polynomial factors of `M_{b+2} = (b−2)(b+2)N₂/6` and `M_{b+1} = (b−2)(a−b+1)N₁/12`.

> **Both standalone bounds are FALSE.** A sweep `m ≤ 45` finds **42 violations each**; e.g.
> `λ=(24,7,3)`, `m=17`: `v₂(N₂)=1` but `v₂(P+1)=3`, so `v₂(N₂)=1 < 2 = v₂(P+1)−1`.

Yet `Δ(b+i) > −θ` never fails. The resolution: the `v₂(N_i)` deficit is **compensated** inside the
descent chain by the surplus the top index carries. So one must *not* isolate `v₂(N_i)`; one bounds
the whole `Δ(b+i)` at once. This is exactly the "deficit is not a single product factor" obstruction
named in PROVE.md — and the fix is to keep the deficit attached to its consecutive-product surplus.

---

## 1. The direct `Δ` formulas

Rather than telescope through the messy `N_i`, we compute each `Δ(b+i)` **directly** from the hook
formula (`f = M_0`), the boundary closed forms `M_{b+i}`, Legendre (`v₂(k!) = k − s₂(k)`), and Kummer
(`carries(x,y) = s₂(x)+s₂(y)−s₂(x+y) = v₂C(x+y,x)`). The factorial linear-parts cancel identically
(using `a − m − P = 0`); only digit-sums and the `N_i`/deficit valuations survive.

> **Lemma D (direct boundary valuations).** For `λ=(a,b,3)` in the box-interior range,
> `Δ(b+3) = (b+5) − 2 s₂(b+5) + 2\,carries(2P,\,b+5) − 2 v₂(a−1) − 2 v₂(P+2)`,
> `Δ(b+2) = (b+2) − 2 s₂(b+3) + 2\,carries(2(P+1),\,b+3) + 2 v₂(N₂) − 2 v₂(a−1) − 2 v₂(P+2)`,
> `Δ(b+1) = (b−1) − 2 s₂(b+1) + 2\,carries(2(P+2),\,b+1) + 2 v₂(N₁) − 2 v₂(a−1)`.

*(Verified against Murnaghan–Nakayama for all `λ=(a,b,3)`, `m ≤ 40`: 0 mismatch. The `N_i` closed
forms themselves cross-checked `m ≤ 29`.)* The first line is the general-`a` Lemma T of the companion
note; the `a`-even specialisation (`v₂(a−1)=0`) recovers its §4.1.

Two useful exact factorisations, obtained by reducing `N_i` modulo `P+1, P+2` (`a = 2P+b+3`):

> `N₂ = 2\big[(P+1)(b+3) + b\big]`,
> `N₁ = (P+1)\cdot 2(b+2)(b+3) + b(5b+7) = (P+2)\cdot 2(b+2)(b+3) + 3(b²−b−4)`.

---

## 2. The factor-in-product principle, one and two factors

> **Lemma F1 (one factor).** Let `R ≥ 0`, `Q ≥ 4`. Among the `Q` consecutive integers
> `R+1,…,R+Q`, for any designated index `t₀`,
> `  v₂\big(\textstyle\prod_{t=1}^Q (R+t)\big) − v₂(R+t₀) \ge 1.`

> **Lemma F2 (two factors).** Let `R ≥ 0`, `Q ≥ 6`. Among the `Q` consecutive integers
> `R+1,…,R+Q`, for any two **distinct** designated indices `t₁ ≠ t₂`,
> `  v₂\big(\textstyle\prod_{t=1}^Q (R+t)\big) − v₂(R+t₁) − v₂(R+t₂) \ge 1.`

*Proof.* The left side equals `v₂` of the product of the remaining factors,
`∏_{t≠t₀}` resp. `∏_{t∉\{t₁,t₂\}}`. Among `Q` consecutive integers exactly `⌊Q/2⌋` or `⌈Q/2⌉` are
even, so at least `⌊Q/2⌋` are. For F1 (`Q ≥ 4`): `⌊Q/2⌋ ≥ 2`, so even after dropping the one factor
`R+t₀` an even factor remains, contributing `≥ 1`. For F2 (`Q ≥ 6`): `⌊Q/2⌋ ≥ 3`, so after dropping
two factors an even factor remains. ∎

The threshold `Q ≥ 6` in F2 is sharp: at `Q = 5` with an odd start there are only two even factors,
and removing both (e.g. `R=0`, `t₁=2`, `t₂=4`) leaves an odd product. *(F1 verified `Q ≤ 39`,
`R < 600`; F2 verified `Q ≤ 39`, `R < 600`, all pairs; `Q=5` failures exhibited — 0/100 spurious.)*

F2 is the **new engine**: it absorbs **two simultaneous unit deficits** `v₂(R+t₁), v₂(R+t₂)` into one
consecutive product. The companion note's single-factor mechanism is the `Q≥4` special case.

---

## 3. `a` even (`θ = 0`): three clean single-factor bounds

Here `b` is odd `≥ 7`, `v₂(a−1) = 0`, and the digit-sum terms collapse via `b+2k` even. We record
the immediate fact:

> **`a` even ⟹ `v₂(N₂) = 1`.** Indeed `b+3` is even and `b` is odd, so `(P+1)(b+3) + b` is **odd**;
> thus `N₂ = 2·(\text{odd})`. ∎

**Top `j = b+3`.** (Companion §4.1.) With `Q' := (b+5)/2 ≥ 6`, `b+5 = 2Q'`,
`carries(2P,b+5) = carries(P,Q')`, and `(b+5) − 2s₂(b+5) = 2v₂(Q'!)`, Lemma D gives
`Δ(b+3) = 2\big[v₂(\prod_{t=1}^{Q'}(P+t)) − v₂(P+2)\big]`. As `P+2` is the `t=2` factor and `Q' ≥ 6`,
F1 yields `Δ(b+3) ≥ 2 > 0`.

**Subtop `j = b+2`.** Put `Q₂ := (b+3)/2 ≥ 5`. Then `b+3 = 2Q₂`,
`carries(2(P+1),b+3) = carries(P+1,Q₂)`, `(b+3)−2s₂(b+3) = 2v₂(Q₂!)`, and `v₂(N₂)=1`, so Lemma D
collapses to
> `Δ(b+2) = 1 + 2\big[\,v₂\big(\textstyle\prod_{t=1}^{Q₂}(P+1+t)\big) − v₂(P+2)\,\big].`
The product is `(P+2)(P+3)\cdots(P+1+Q₂)`; `P+2` is its `t=1` factor, so the bracket is `≥ 0` (it is
just "a factor divides the product"), giving `Δ(b+2) ≥ 1 > 0`.

**Subtop `j = b+1`.** Put `Q₁ := (b+1)/2 ≥ 4`. Then `b+1 = 2Q₁`,
`carries(2(P+2),b+1) = carries(P+2,Q₁)`, `(b+1)−2s₂(b+1) = 2v₂(Q₁!)`, so
> `Δ(b+1) = 2\big[\,v₂\big(\textstyle\prod_{t=1}^{Q₁}(P+2+t)\big) + v₂(N₁) − 1\,\big].`
The product is `(P+3)(P+4)\cdots(P+2+Q₁)`, **`Q₁ ≥ 4` consecutive integers**, among which two are
even and one is `≡ 0 \pmod 4`; hence its valuation is `≥ 3`. With `v₂(N₁) ≥ 0` the bracket is
`≥ 3 + 0 − 1 = 2`, so `Δ(b+1) ≥ 4 > 0`.

So **all three `a`-even boundary indices satisfy `Δ > 0`**, with margins `≥ 2, 1, 4` respectively.

---

## 4. `a` odd (`θ = 3`): three bounds via the two-factor engine

Now `b = 2β` is even, `β ≥ 5`. The two deficits `v₂(a−1) = 1 + v₂(P+β+1)` and
`v₂(P+2)` (from `a−b+1`) are **both positive** in general — the genuinely harder regime. The carries
collapse uses the odd offsets (`b+2k+1` odd): writing `b+5 = 2(β+2)+1` etc., the bit-0 of the odd
summand never carries, so
`carries(2P, b+5) = carries(P, β+2)`, `carries(2(P+1), b+3) = carries(P+1, β+1)`,
`carries(2(P+2), b+1) = carries(P+2, β)`, and `(b+2k+1) − 2s₂(b+2k+1) = 2v₂((β+\cdots)!) − (\text{odd})`.

Substituting into Lemma D and folding factorial-valuation `+` carries into one consecutive product via
Kummer (`v₂(j!) + carries(R,j) = v₂(\prod_{t=1}^{j}(R+t))`), and using
`v₂(a−1) = 1 + v₂(P+β+1)`, `v₂(N₂) = 1 + v₂(n₂)` with `n₂ := (P+1)(2β+3) + 2β`:

> `Δ(b+3) = 2\big[\,U   − v₂(P+β+1) − v₂(P+2)\,\big] − 3`, `U   := v₂\big(\prod_{t=1}^{β+2}(P+t)\big) = v₂\big((P+1)\cdots(P+β+2)\big)`,
> `Δ(b+2) = 2\big[\,U'  − 1 + v₂(N₂) − v₂(a−1) − v₂(P+2)\,\big]`, `U'  := v₂\big((P+2)\cdots(P+β+2)\big)`,
> `Δ(b+1) = 2\big[\,U'' + v₂(N₁) − v₂(a−1)\,\big] − 3`, `U'' := v₂\big((P+3)\cdots(P+β+2)\big)`.

*(All three identities verified against Murnaghan–Nakayama for all `a`-odd `λ=(a,b,3)`, `m ≤ 41`:
0 mismatch.)* The three consecutive products have lengths `β+2 ≥ 7`, `β+1 ≥ 6`, `β ≥ 5` respectively,
and each contains the relevant deficit factors:

- `U` contains `P+2` (`t=2`) and `P+β+1` (`t=β+1`), **distinct** since `β ≥ 5`;
- `U'` contains `P+2` (`t=1`) and `P+β+1`, distinct;
- `U''` contains `P+β+1` (`t=β−1`); it does **not** contain `P+2` — and indeed the `Δ(b+1)` formula
  carries only the single deficit `v₂(a−1)` (the `a−b+1` deficit cancelled in the descent).

**Top `j = b+3`** (this is claim (C)). By **F2** on the `β+2 ≥ 7 ≥ 6` consecutive factors of `U`,
removing the two distinct factors `P+2, P+β+1`:
`U − v₂(P+β+1) − v₂(P+2) ≥ 1`. Hence `Δ(b+3) ≥ 2·1 − 3 = −1 > −3`. ✓

**Subtop `j = b+2`.** By **F2** on the `β+1 ≥ 6` consecutive factors of `U'`:
`U' − v₂(P+β+1) − v₂(P+2) ≥ 1`. Therefore the bracket is
`(U' − v₂(P+β+1) − v₂(P+2)) − 1 + v₂(N₂) ≥ 1 − 1 + v₂(N₂) = v₂(N₂) ≥ 1`, so actually
`Δ(b+2) ≥ 2 ≥ 0 > −3`. ✓ *(Using `v₂(N₂) = 1 + v₂(n₂) ≥ 1`; even `v₂(N₂) ≥ 0` already gives
`Δ(b+2) ≥ 0`.)*

**Subtop `j = b+1`.** First, **`a` odd ⟹ `v₂(N₁) ≥ 1`**: in `N₁ = (P+1)·2(b+2)(b+3) + b(5b+7)` with
`b` even, the first summand has `v₂ ≥ 2` (the factor `2(b+2)` with `b+2` even) and the second has
`v₂(b(5b+7)) = v₂(b) ≥ 1` (`5b+7` odd), so `N₁` is even. Next, by **F1** on the `β ≥ 5 ≥ 4`
consecutive factors of `U''`, dropping the single factor `P+β+1`:
`U'' − v₂(P+β+1) ≥ 1`. Hence the bracket is
`(U'' − v₂(P+β+1)) + v₂(N₁) − 1 ≥ 1 + 1 − 1 = 1` (here `v₂(a−1) = 1 + v₂(P+β+1)`), giving
`Δ(b+1) ≥ 2·1 − 3 = −1 > −3`. ✓

So **all three `a`-odd boundary indices satisfy `Δ(b+i) > −3`** (margins: `−1, 0, −1` against the
floor `−3`, i.e. strict by `≥ 2`).

---

## 5. Conclusion

> **Theorem (three-row `c = 3` boundary lemma).** For `λ=(a,b,3) ⊢ 2m` with `a+b = 2m−3`, in the
> box-interior range (`b ≥ 7` for `a` even, `b ≥ 10` for `a` odd), every boundary index
> `j ∈ {b+1,b+2,b+3}` satisfies `val(j) > V`. Equivalently, no boundary index is an extra minimizer
> of `val`, so `J*` equals the interior 2-adic box.

Combined with the interior theorem (`2026-06-14`, `2026-06-15`), the smaller-`b` boundary-tie shapes
handled in the interior note, and the `a`-even-top already in the companion note:

> **Three-row `d=4` family `(a,b,3)` is a COMPLETE theorem.** `J*` is exactly the interior
> 2-adic box `j₀ + \{0,2,4,6\}` (`j₀ ∈ \{0,3\}`); on every tie `|J*| ∈ \{1,2,4\}` is even, so the
> leading-`π` dichotomy gives `G_λ(i) = 0` only on `|J*|`-even shapes — matching `c=1, c=2`. No
> machine residual remains for `c ≤ 3`.

The **two-factor Lemma F2** is the keystone, and it is uniform in `c`: the `a`-odd top index always
presents two simultaneous unit deficits (`v₂(a−c+2)` and `v₂(a−b+1)` both positive), exactly the
configuration F2 dissolves. This is the first instance of the general-`c` `a`-odd-top obstruction and
its general resolution.

### Bonus — the general-`c` boundary lemma (stated)

The same mechanism gives the shape of the general-`c` proof:

> **General-`c` boundary conjecture (now a theorem template).** For `λ=(a,b,c)`, each boundary
> `Δ(b+i)` reduces by Lemma D's recipe to `2[\,v₂(\text{a consecutive product of length} \sim
> (b+\cdots)/2) − (\text{deficits}) + v₂(N_i^{(c)})\,] − (\text{small odd})`. The deficits are exactly
> the factors `P+t₀` of that product coming from the linear forms `a−c+2, a−b+1, P+1, \ldots, P+c`;
> the `a`-even top needs **F1**, the `a`-odd top needs **F2** (two deficits), and the lower indices
> inherit a spare `v₂(N_i^{(c)}) ≥ 0`. The only `c`-specific input is the parity bookkeeping of
> `a−c+2` and the `N_i^{(c)}` (which are even of controlled valuation). F1 + F2 are `c`-uniform.

---

## 6. Verification summary

All scripts in `projects/code/threerow-boundary/` (shared `mn.py` = Murnaghan–Nakayama).

| claim | range | result |
|---|---|---|
| `M_{b+1},M_{b+2},M_{b+3}` closed forms (incl. `N₁,N₂`) vs MN | `m ≤ 29` | 0 mismatch |
| Lemma D (three direct `Δ` formulas) vs MN | `m ≤ 40` | 0 mismatch |
| standalone bounds (A),(B) — **shown FALSE** | `m ≤ 45` | 42 violations each |
| `a`-even product forms `Δ(b+1),Δ(b+2)`; `v₂(N₂)=1`, `W₁≥3` | `m ≤ 41` | 0 mismatch / asserts pass |
| `a`-odd product forms `Δ(b+1),Δ(b+2),Δ(b+3)` | `m ≤ 41` | 0 mismatch |
| `a`-odd key inequalities (F2 top/subtop, F1 subtop, `v₂(N₁)≥1`) | `m ≤ 41` | 0 failure |
| Lemma F1 (`Q≥4`), F2 (`Q≥6`); `Q=5` sharpness | `Q≤39`, `R<600` | 0 failure / `Q=5` fails |
| **end-to-end** `Δ(b+i) > −θ`, both parities, box-interior | `m ≤ 60` (3828 indices) | **0 violation** |
| boundary scout strict inequality (prior CODE session) | `m ≤ 200` (~59k) | 0 violation |

### Files
`projects/code/threerow-boundary/`: `explore_N.py`, `explore_N1.py` (the `N_i` factorisations),
`direct_forms.py` (Lemma D), `aeven_check.py` (§3 forms), `aodd_check2.py` (§4 forms + key
inequalities), `lemmaF2.py` (Lemmas F1/F2 + sharpness), `mn_crosscheck.py` (`M` closed forms),
`final_endtoend.py` (boundary lemma, `m ≤ 60`).

### Gaps
None for `c ≤ 3`. The general-`c` template (§5 bonus) is stated but not written out for `c ≥ 4`;
that is the natural next target (and the LEAN target alongside the now-complete `c=1,2,3` families).
