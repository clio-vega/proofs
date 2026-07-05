# The `c = 6` interior is a genuine two-generator even-|J*| theorem: the mod-4 table, proved

**Date:** 2026-07-04 (prove session, part II)
**Status:** **PROVED (no residual).** The box-interior three-row family `λ = (a,b,6)` has minimiser
set `J*` a deterministic function of `(a,b) mod 4`, with `|J*| ≤ 2` even on every tie:

> **Theorem (c = 6 interior).** For every box-interior three-row shape `λ = (a,b,6)`
> (`a = 2P + b + 6`, `b ≥ 12`, so `a ≡ b (mod 2)`), the interior minimiser set is
>
> | `(a,b) mod 4` | `J*` | firing generator |
> |---|---|---|
> | `(0,0)`, `(3,1)` | `{0, 2}` | 2 |
> | `(0,2)`, `(1,1)` | `{0, 4}` | **4** |
> | `(1,3)`, `(2,0)`, `(2,2)`, `(3,3)` | `{0}` | — |
>
> In every case `0 ∈ J*`, `J* ⊆ {0, 2, 4}`, and `|J*| ≤ 2` is **even**. Hence
> `G_λ(i) = 0 ⟺ |J*|` even holds with no residual for the `c = 6` interior.

This makes `c = 6` the **first even `c > 2` proved as a complete interior theorem with two genuine
generators** (generators 2 and 4 both fire, on disjoint `(a,b) mod 4` classes). It confirms
Conjecture 1′ of `2026-07-04-c6-generator4-refutes-single-generator.md` (which refuted the
single-generator law by exhibiting `J* = {0,4}` at `(13,13,6)`) and turns it into a proof.

The crux — *why generator 4 fires for `c = 6` but not for `c = 4`* — has a clean and, it turns out,
**surprising** answer: it is an arithmetic resonance in `c mod 4`, **not** the run-length effect
speculated in the refutation note. §5.

All notation is that of `2026-06-19-c4-interior-number-lemma.md` and
`2026-06-18-c4-interior-Jstar-even.md`, whose engine transfers to `c = 6` verbatim.

---

## 0. Setup, closed form, decomposition, Prop 2

With `h₁ = p₁`, `e₂ = s_{(1,1)}`, `M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩`, `val(j) = j + 2 v₂(C(m,j) M_j)`,
`J* = {j : val(j) = min}`, and `Δ(j) := val(j) − val(0)`. Leading-`π` dichotomy (proved earlier,
Lean-checked): `C ≡ |J*| (mod π)`, so `|J*|` odd `⟹ G_λ(i) ≠ 0`; ties are exactly `|J*|` even.
Since `a + b = 2m − 6` is **even**, `a ≡ b (mod 2)`, so `j₀ = 0` and — because `val(j) ≡ j (mod 2)`
— **only even `j` can tie with `j = 0`.**

Write `s₂` for the binary digit sum, `vfact(r) = v₂(r!) = r − s₂(r)`, `(x)_r = x(x−1)⋯(x−r+1)`.

**Closed form** (`c`-uniform, verified `c = 6` vs Murnaghan–Nakayama, `c6_build.py`, 6349 shapes,
0 mismatch):
```
M_j = C(2(m−j), b−j) (a−b+1) Q_6 / [ 6! (a+7−j) ∏_{i=1}^{6}(b+i−j) ],
Q_6(a,b,j) = (a−4)(b−5) H_6(a,b,j) + 12! · C(j,12),
```
with the **inhomogeneous tip** `12!·C(j,12) = (j)_{12}` vanishing for `0 ≤ j ≤ 11`, the **heavy
factor** `L_6 = (a−4)(b−5)` (the `c = 6` instance of `(a−(c−2))(b−(c−1))`), and `H_6` of degree
`10` in `j`. Its anchor value is the product of two length-`5` consecutive runs:
```
H_6(0) = (a+3)(a+4)(a+5)(a+6)(a+7) · (b+2)(b+3)(b+4)(b+5)(b+6).
```

**Proposition 2 (c = 6).** For `1 ≤ j ≤ b`,
```
Δ(j) = j − 2 s₂(j) + 2 v₂C(a+7, j) + 2 v₂C(b+6, j) + 2[ v₂Q_6(j) − v₂Q_6(0) ].
```
(Verified vs direct `val(j) − val(0)` over MN, `c6_residue.py`, 3079 cases, 0 mismatch.) Write
`T(j) = v₂C(a+7,j) + v₂C(b+6,j) + v₂Q_6(j) − v₂Q_6(0)`, so `Δ(j) = j − 2 s₂(j) + 2 T(j)`.

**Peeling identity.** The top `c−1 = 5` factors of `(a+7)_j` are `(a+7)(a+6)(a+5)(a+4)(a+3)` — the
`a`-run of `H_6(0)` — and the top 5 of `(b+6)_j` are `(b+6)(b+5)(b+4)(b+3)(b+2)` — the `b`-run. They
cancel `H_6(0)`. For `j ≥ 5`, with the remaining runs `(a+2)_{j−5}, (b+1)_{j−5}` (verified
`verify_peeling_caseA.py`):
```
T(j) = v₂((a+2)_{j−5}) + v₂((b+1)_{j−5}) − 2 vfact(j) + v₂Q_6(j) − v₂(a−4) − v₂(b−5).
```

**The Number Lemma N6.** `content(H_6) = 7` over the parity sheet `a ≡ b (mod 2)`, attained at
**every** `j`-slice (per-`j` content `= 7` for `j = 0,…,12,16`, `c6_build.py`). Equivalently:

> **Lemma N6.** `2⁷ ∣ H_6(a,b,j)` for all `a ≡ b (mod 2)` and all `j` (`v₂H_6(j) ≥ 7`), and this
> constant floor `β'(6) = 7` is sharp.

*Proof.* `H_6` is an integer polynomial, so `H_6 mod 2⁷` depends only on `(a,b,j) mod 2⁷`; the
2-adic content over the parity sheet is exactly the minimum 2-adic valuation of the coefficients in
the integer-valued (forward-difference) basis, computed exactly in `c6_build.py` / `jobB` —
`= 7`. ∎ (This is the `c = 6` analogue of the `c = 4` Number Lemma `16 ∣ H`, and equals the anchor
content `γ(6) = 7` of the run-content theorem, `2026-07-04-...refutes...`, Theorem A, since the dip
`γ(6) − β'(6) = 0`.)

For `1 ≤ j ≤ 11` the tip vanishes, so `Q_6(j) = (a−4)(b−5) H_6(j)`, and `v₂Q_6(j) − v₂Q_6(0)`
`= v₂H_6(j) − v₂H_6(0)`. In particular the peeling identity collapses on `5 ≤ j ≤ 11` to
```
T(j) = v₂((a+2)_{j−5}) + v₂((b+1)_{j−5}) − 2 vfact(j) + v₂H_6(j).            (★)
```

The proof proceeds by régime: low `j ∈ {1,2,3,4}` (below peeling — the two generators live here),
mid `j ∈ {5,…,11}` (peeling, tip-free), deep `j ≥ 12` (tip active).

---

## 1. Low indices `j = 1, 2, 3, 4` — the two generators

For `1 ≤ j ≤ 4`, cancelling the shortened runs against the top factors of `C(a+7,j), C(b+6,j)`
(a bookkeeping identical to `c = 4` §3, verified symbolically) collapses `T(j)` to `v₂` of a single
**leftover polynomial** `G_j := H_6(j) / (\text{shortened runs})`. Writing out the cancellation for
each `j` gives (constants are `−2·[v₂ of the binomial denominators cancelled]`):

| `j` | `Δ(j)` in terms of `G_j` | requirement for `Δ(j) ≥ 0` |
|----|----|----|
| 1 | `Δ(1) = 2 v₂(G_1) − 1` | `v₂(G_1) ≥ 1` |
| 2 | `Δ(2) = 2 v₂(G_2) − 4` | `v₂(G_2) ≥ 2` |
| 3 | `Δ(3) = 2 v₂(G_3) − 5` | `v₂(G_3) ≥ 3` |
| 4 | `Δ(4) = 2 v₂(G_4) − 10` | `v₂(G_4) ≥ 5` |

with
```
G_1 = a b + 6a + 7b + 12 = (a+7)(b+6) − 30,
G_2 = a²b² + 11a²b + 30a² + 13ab² + 71ab + 78a + 42b² + 78b + 36,
G_3 = H_6(3)/[(a+3)(a+4)(b+2)(b+3)]        (bicubic),
G_4 = H_6(4)/[(a+3)(b+2)]                   (biquartic).
```

**`j = 1` (proved, `Δ > 0`).** `G_1 = (a+7)(b+6) − 30`. Since `a ≡ b (mod 2)`, exactly one of
`a+7, b+6` is even, so `(a+7)(b+6)` is even; with `30` even, `G_1` is even, `v₂(G_1) ≥ 1`, and
`Δ(1) ≥ 1 > 0`. ∎ (`j = 1` is odd, so `val(1)` is odd and cannot tie; we need only `Δ(1) > 0`.)

**`j = 3` (proved, `Δ > 0`).** The exact 2-adic content of `G_3` over the parity sheet is `3`
(`c6_residue.py`, binomial-basis, rigorous), i.e. `2³ ∣ G_3` always. Hence `Δ(3) = 2 v₂(G_3) − 5
≥ 1 > 0`. ∎ (`j = 3` odd, cannot tie.)

**`j = 2` — generator 2 (proved).** The exact content of `G_2` over the parity sheet is `2`, so
`2² ∣ G_2` always and `Δ(2) = 2 v₂(G_2) − 4 ≥ 0`: `0 ∈ J*` at level 2. Computing the content of
`G_2` restricted to each `(a,b) mod 4` class (exact binomial basis on the sublattice `4ℤ+r`,
`c6_certify.py`) gives
```
content(G_2) = 2  on classes (0,0), (3,1);      content(G_2) = 3  on the other 6 classes.
```
On the six classes with content `3`, `Δ(2) ≥ 2·3 − 4 = 2 > 0`. On the two classes `(0,0),(3,1)`,
`v₂(G_2)` is **constant equal to 2** (exhaustive mod-8 check `c6_tieconst.py`: over every mod-8
lift of these two classes `v₂(G_2) = 2`), so `Δ(2) = 0` there. Therefore
```
Δ(2) = 0  ⟺  (a,b) ≡ (0,0) or (3,1)  (mod 4),      else  Δ(2) ≥ 2.     ∎
```

**`j = 4` — generator 4 (proved).** The exact content of `G_4` over the parity sheet is `5`, so
`2⁵ ∣ G_4` always and `Δ(4) = 2 v₂(G_4) − 10 ≥ 0`: `0 ∈ J*` at level 4 too. Per-class content
(`c6_certify.py`):
```
content(G_4) = 5  on classes (0,2), (1,1);      content(G_4) = 6  on the other 6 classes.
```
On the six content-`6` classes, `Δ(4) ≥ 2·6 − 10 = 2 > 0`. On `(0,2),(1,1)`, `v₂(G_4)` is
**constant equal to 5** (exhaustive mod-64 check `c6_tieconst.py`), so `Δ(4) = 0` there. Therefore
```
Δ(4) = 0  ⟺  (a,b) ≡ (0,2) or (1,1)  (mod 4),      else  Δ(4) ≥ 2.     ∎
```

**Disjointness.** The generator-2 classes `{(0,0),(3,1)}` and generator-4 classes `{(0,2),(1,1)}`
are **disjoint**. So at most one of `Δ(2), Δ(4)` vanishes, giving `|J* ∩ {0,2,4}| ≤ 2` and — once
§§2–4 show every other `Δ(j) > 0` — exactly the mod-4 table. Because `Δ(2), Δ(4) ≥ 0` always,
`min_j val(j) = val(0)`, so **`0 ∈ J*` unconditionally.**

---

## 2. Mid indices `j = 5, …, 11` (peeling, tip-free) — proved `Δ(j) > 0`

Here `(★)` holds: `Δ(j) = j − 2 s₂(j) + 2[ v₂(Φ_j) − 2 vfact(j) ]` with
`Φ_j := (a+2)_{j−5} (b+1)_{j−5} H_6(j)`. Each `Δ(j) > 0` is thus the finite divisibility
`v₂(Φ_j) ≥ 2 vfact(j) − (j − 2 s₂(j))/2 + 1`, a complete residue-class statement decided by the
exact 2-adic content of the polynomial `Φ_j` over the parity sheet (`c6_residue.py`, binomial
basis — rigorous, no `m`-cap):

| `j` | `content(Φ_j)` | `vfact(j)` | `min Δ(j)` |
|----|----|----|----|
| 5 | 7 | 3 | **3** |
| 6 | 8 | 4 | **2** |
| 7 | 9 | 4 | **3** |
| 8 | 12 | 7 | **2** |
| 9 | 14 | 7 | **5** |
| 10 | 15 | 8 | **4** |
| 11 | 16 | 8 | **5** |

All `min Δ(j) ≥ 2 > 0`. ∎ (These are the `c = 6` analogues of the `c = 4` §4 residue checks; the
even indices `6, 8, 10` are the would-be generators — all strictly positive.)

---

## 3. Deep indices `j ≥ 12` (tip active) — proved `Δ(j) > 0`

Fix `j ≥ 12`. Split on the two summands of `Q_6(j) = (a−4)(b−5) H_6(j) + (j)_{12}`. Write
`P := v₂((j)_{12}) = vfact(j) − vfact(j−12) = 12 + s₂(j−12) − s₂(j)` (`≥ 10`), and
`Hv := v₂((a−4)(b−5) H_6(j)) = φ + v₂H_6(j)`, `φ := v₂((a−4)(b−5)) ≥ 1`.

**Case A: `Hv ≥ P`** (so `v₂Q_6(j) ≥ P`). For `j ≥ 12` the run `(a+2)_{j−5}` (length `≥ 7`) splits
as a top 6-block `(a+2)(a+1)a(a−1)(a−2)(a−3)` (`v₂ ≥ vfact(6) = 4`), the factor `a−4`
(`= L_6` on the `a`-side), and a tail `(a−5)_{j−12}` (`v₂ ≥ vfact(j−12)`); symmetrically for `b`.
Substituting in the peeling identity (the `v₂(a−4), v₂(b−5)` cancel):
```
T(j) ≥ 8 + 2 vfact(j−12) − 2 vfact(j) + v₂Q_6(j) = 8 − P + v₂Q_6(j) ≥ 8 − P,
Δ(j) ≥ j − 2 s₂(j) + 16 − 2P = j − 8 − 2 s₂(j−12)  (using P = 12 + s₂(j−12) − s₂(j)).
```
Writing `j = 12 + t`: `Δ(j) ≥ 4 + (t − 2 s₂(t)) ≥ 3`, since `t − 2 s₂(t) ≥ −1` for all `t ≥ 0`
(Lemma A). **So `Δ(j) ≥ 3 > 0`.** ∎ (Case A) — the `c`-uniform bound `2c − 4 s₂(c) − 1 = 3`.

**Case B: `Hv < P`** (tip strictly larger valuation, so `v₂Q_6(j) = Hv`; the peeling identity's
`−φ` cancels the `+φ`). Then `Δ(j) = Δ̂(j)` with
```
Δ̂(j) = j − 2 s₂(j) + 2[ v₂((a+2)_{j−5}) + v₂((b+1)_{j−5}) − 2 vfact(j) + v₂H_6(j) ].
```
Using `v₂((x)_{j−5}) ≥ vfact(j−5)` and Lemma N6 `v₂H_6(j) ≥ 7`:
```
Δ̂(j) ≥ g_6(j) := j − 2 s₂(j) − 4 v₂((j)_5) + 14 = j + 2 s₂(j) − 4 s₂(j−5) − 6,
```
(the last equality via `v₂((j)_5) = 5 − s₂(j) + s₂(j−5)`).

*Claim: `g_6(j) > 0` for all `j ≥ 12` except `j ∈ {12, 16}`.* For `j ≥ 32`,
`g_6(j) ≥ j − 8 − 4 log₂ j > 0` (drop `2 s₂(j) ≥ 2`, bound `s₂(j−5) ≤ log₂ j + 1`; the real
function `x − 8 − 4 log₂ x` is increasing for `x > 4/\ln 2` and `= 4` at `x = 32`). For
`12 ≤ j ≤ 31` the values are computed directly (`c6_certify.py`): `g_6(j) > 0` **except**
`g_6(12) = −2` and `g_6(16) = 0`. ∎ (claim)

**The two deep exceptions `j ∈ {12, 16}`.** Here `Δ̂(j) = j − 2 s₂(j) + 2[ v₂(Π_j) − 2 vfact(j) ]`
with `Π_j := (a+2)_{j−5} (b+1)_{j−5} H_6(j)`. Both `j` even, so `Δ̂(j) > 0 ⟺ Δ̂(j) ≥ 2 ⟺`:
```
j = 12:  Δ̂(12) = 2 v₂(Π_{12}) − 32 ≥ 2  ⟺  2^{17} ∣ Π_{12},
j = 16:  Δ̂(16) = 2 v₂(Π_{16}) − 46 ≥ 2  ⟺  2^{24} ∣ Π_{16}.
```
Each is a polynomial divisibility over `a ≡ b (mod 2)`, decided by the **exact 2-adic content** of
`Π_j` (binomial basis — a complete residue-class certificate, `c6_lowj.py`):
```
content(Π_{12}) = 19 ≥ 17   ⟹  Δ̂(12) ≥ 6 > 0,
content(Π_{16}) = 27 ≥ 24   ⟹  Δ̂(16) ≥ 8 > 0.
```
(The crude independent floors — `vfact(7)+vfact(7)+7 = 15`, `vfact(11)+vfact(11)+7 = 21` — fall
short; the extra 2-content comes from the *joint* correlation of the two runs with `H_6`, exactly as
at `c = 4` `j ∈ {8,10}`. Computing `content(Π_j)` directly captures it without brute force over
`2^{17}` residues.) ∎ (Case B)

Combining: `Δ(j) > 0` for every deep `j ≥ 12`.

---

## 4. Assembly

For every interior `1 ≤ j ≤ b`:
- `Δ(j) > 0` for `j ∈ {1, 3}` (§1), `j ∈ {5,…,11}` (§2), and `j ≥ 12` (§3);
- `Δ(2) ≥ 0`, with `Δ(2) = 0 ⟺ (a,b) ≡ (0,0),(3,1) (mod 4)` (§1);
- `Δ(4) ≥ 0`, with `Δ(4) = 0 ⟺ (a,b) ≡ (0,2),(1,1) (mod 4)` (§1);
- the generator-2 and generator-4 classes are disjoint.

Hence `min_j val(j) = val(0)`, `0 ∈ J*`, `J* ⊆ {0,2,4}`, and `J*` is exactly the mod-4 table.
Every tie has `|J*| = 2` (even), so `G_λ(i) = 0 ⟺ |J*|` even holds with **no residual** for the
`c = 6` interior. ∎

**Cross-validation.** End-to-end Murnaghan–Nakayama census of `J*` for all interior `(a,b,6)`,
`12 ≤ b`, `a ≤ 120` (`c6_residue.py`): **1488 shapes, 0 mismatch** with the mod-4 table.

---

## 5. The crux: why generator 4 fires for `c = 6` — and a correction

Generator 4 ties iff `Δ(4) = 0 ⟺ v₂(G_4^{(c)}) = 5`, where `G_4^{(c)}` is the reduced heavy at
`j = 4` (the leftover after cancelling the shortened runs; `Δ(4) = 2 v₂(G_4^{(c)}) − 10` for every
even `c ≥ 6`, the constant `10` coming from `4 − 2s₂(4) − 4·v₂(4!) = 2 − 12`). The tie is thus a
statement about **one unit** of 2-adic content: `content(G_4^{(c)})` equals `5` (tie) or `6`
(strict). Direct Murnaghan–Nakayama computation of `min Δ(4)` (`c6_crux.py`) gives:

| `c` | 4 | 6 | 8 | 10 | 12 |
|----|----|----|----|----|----|
| `min Δ(2)` (gen 2) | 0 | 0 | 0 | 0 | 0 |
| `min Δ(4)` (gen 4) | **2** | **0** | **2** | **0** | **2** |

So **generator 4 fires iff `c ≡ 2 (mod 4)`** (`c = 6, 10`), and is suppressed for `c ≡ 0 (mod 4)`
(`c = 4, 8, 12`); generator 2 fires for every even `c`. This is the precise resolution of the
`c=4`-vs-`c=6` puzzle.

**A correction to the refutation note.** `2026-07-04-c6-generator4-refutes-single-generator.md` §4
attributed the appearance of generator 4 to the *run length* `c − 1` ("the anchor's run length sets
the size of the generator menu … generators above 2 proliferate as `c` grows"). The data above
**refutes** the monotone reading: the menu does **not** grow with `c`. `c = 8` (runs of length 7,
longer than `c = 6`'s length 5) has **no** generator 4. The true governor is the parity `c mod 4`
of the reduced-heavy content at `j = 4`, an arithmetic resonance, not a size effect. (What survives
of §4 is only that the tie is a low-régime, tip-free `H_c`-content event — correct — but its
*incidence* is periodic in `c`, not monotone.)

For `c = 4`, `j = 4` sits *at* the peeling threshold `c − 1 = 3` (`4 > 3`), so its analysis is the
`c = 4` §4 residue check `2⁶ ∣ (a+2)(b+1)H_4(4)`, whose content is `6` (min `Δ(4) = 2`) — no tie.
For `c = 6`, `j = 4` is *below* threshold `c − 1 = 5`, a pure cancellation leaving `G_4` of content
`5` on `(0,2),(1,1)` — the tie. The `c ≡ 2 mod 4` law predicts (and MN confirms `c ≤ 12`) that this
one-unit content dip recurs with period 4 in `c`; a `c`-uniform proof of that periodicity is the
natural next question (it would classify generator 4 for all even `c` at one stroke).

---

## 6. What is proved

**Proved in full (symbolic / exact binomial-basis content — complete residue-class certificates,
no `m`-cap):**
- Closed form + decomposition `Q_6 = (a−4)(b−5)H_6 + (j)_{12}`; Prop 2; the peeling identity
  (verified vs MN, 0 mismatch throughout).
- **Lemma N6:** `2⁷ ∣ H_6` (constant floor `β'(6) = 7`, per-`j` content `= 7`, sharp).
- Low `j`: `Δ(1),Δ(3) > 0`; `Δ(2) ≥ 0` with exact mod-4 tie class `{(0,0),(3,1)}`; `Δ(4) ≥ 0` with
  exact mod-4 tie class `{(0,2),(1,1)}` (contents `2,3` resp. `5,6`; tie-class constancy by
  exhaustive mod-8 / mod-64 checks).
- Mid `j = 5,…,11`: `Δ(j) ≥ 2 > 0` (exact content of `Φ_j`).
- Deep `j ≥ 12`: Case A `Δ ≥ 3`; Case B `Δ̂ ≥ g_6(j) > 0` off `{12,16}` (analytic `j ≥ 32`, direct
  `12 ≤ j ≤ 31`); `j ∈ {12,16}` closed by `content(Π_{12}) = 19 ≥ 17`, `content(Π_{16}) = 27 ≥ 24`.

**Established (data, `c ≤ 12`, clearly labelled as beyond this theorem):**
- Generator 4 fires iff `c ≡ 2 (mod 4)` (min `Δ(4)`: `2,0,2,0,2` for `c = 4,6,8,10,12`), correcting
  the run-length heuristic of the refutation note.

**Net.** The `c = 6` interior even-|J*| law is a **complete theorem**: `J*` is the mod-4 table,
`|J*| ≤ 2` even, no residual — the first even `c > 2` with two genuine generators. Its interior now
stands with `c = 1,2,3,4,5`; only the `c = 6` boundary remains for full family completion (the
general-`c` boundary master of `2026-06-17-...` supplies `θ(6) = 0` and the content route).

### For LEAN
`2⁷ ∣ H_6` is `decide`-friendly exactly as `16 ∣ H` was: a finite residue check mod `2⁷` on a
polynomial with small integer coefficients. The whole interior reduces to a fixed finite list of
divisibilities `2^{k_j} ∣ Φ_j` (for `j ≤ 16`) plus the analytic `g_6(j) > 0` tail and the `c`-uniform
Case A / Lemma A arithmetic — all decidable or already-Lean lemmas (`vz_choose_ge`, Kummer).

### Files (`projects/code/threerow-heavy/`)
`c6_build.py` (H_6 build, decomposition, per-`j` content), `c6_lowj.py` (leftover factors + deep
exceptions `Π_{12},Π_{16}`), `c6_residue.py` (Prop-2 vs MN, low-`j`/mid-`j` contents, MN census),
`c6_certify.py` (exact per-mod-4-class contents + `g_6` tail), `c6_tieconst.py` (mod-8 / mod-64
tie-class constancy), `c6_crux.py` (`min Δ(4)` across `c = 4..12`). MN engine `threerow-c3/mn.py`;
`H_c` reconstruction `jobB_betaprime_sequence.py`.
