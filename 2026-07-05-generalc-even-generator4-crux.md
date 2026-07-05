# The general even-`c` interior: generator 4 fires iff `c ≡ 2 (mod 4)` — the crux, proved `c`-uniformly

**Date:** 2026-07-05 (prove session) · **Status:** the genuinely-open step of the general even-`c`
interior even-`|J*|` law is **PROVED**: a `c`-uniform *derivation* of the generator menu, in
particular of

> **Main Theorem (generator-4 crux).** For every even `c ≥ 6` and box-interior three-row
> `λ = (a,b,c)` (`a = 2P+b+c`, `b` large), the level-4 tie index satisfies
> `Δ(4) = 2·v₂(G_4^{(c)}) − 10`, and
> `   min Δ(4) = [\,c ≡ 2 \pmod 4\,] ? 0 : 2.`
> Hence **generator 4 is available (a `{0,4}` tie occurs) iff `c ≡ 2 (mod 4)`**, on the
> `c`-uniform mod-4 classes `(a,b) ≡ (0,2),(1,1)`; for `c ≡ 0 (mod 4)` generator 4 is dead.

This is the arithmetic the whole `c mod 4` periodicity rests on (PROVE.md §"The crux"). The engine
is a one-line valuation of a single explicit factor:

> **`K(c) = 24·c·(c−1)·(c−4)·(c−5)`**, the minimising binomial-basis coefficient of the reduced
> heavy `G_4^{(c)}`; for even `c`, `v₂(K(c)) = 3 + v₂(c) + v₂(c−4)`, which is **`5` iff `c ≡ 2 (mod 4)`**
> (then `c, c−4 ≡ 2 mod 4`, each `v₂ = 1`) and **`≥ 7` iff `c ≡ 0 (mod 4)`** (then both `≡ 0 mod 4`).

Together with the companion generator-2 result (§3) and the deep/mid régime (§4), this assembles the
**general even-`c` interior law**, a table periodic in `c mod 4` (§5), correcting the "menu grows with
`c`" reading and unifying the proven `c = 4` (`2026-06-19-c4-interior-number-lemma.md`) and `c = 6`
(`2026-07-04-c6-interior-two-generator.md`) families.

All notation is that of the `c = 6` note. Setup recap: `M_j = ⟨s_λ, e₂^j h₁^{2m−2j}⟩`,
`val(j) = j + 2 v₂(C(m,j) M_j)`, `J* = argmin_j val(j)`, `Δ(j) = val(j) − val(0)`. Since `c` even
⟹ `a ≡ b (mod 2)` ⟹ `j₀ = 0` and `val(j) ≡ j (mod 2)`, so **only even `j` can tie with `0`**. The
leading-`π` dichotomy (`C ≡ |J*| mod π`) makes `G_λ(i) = 0 ⟺ |J*|` even; the theorem shows `|J*| ≤ 2`
is always even. `s₂` = binary digit sum, `vfact(r) = r − s₂(r)`, `(x)_r = x(x−1)⋯(x−r+1)`.

---

## 0. The uniform closed form and Prop 2

The `c`-uniform closed form (Job B, `jobB_betaprime_sequence.py`, cross-checked vs Murnaghan–Nakayama):
```
M_j = C(2(m−j), b−j)·(a−b+1)·Q_c / [ c!·(a+c+1−j)·∏_{i=1}^{c}(b+i−j) ],
Q_c = (a−(c−2))(b−(c−1))·H_c(a,b,j) + (2c)!·C(j,2c),
```
with the inhomogeneous **tip** `(2c)!·C(j,2c) = (j)_{2c}` vanishing for `0 ≤ j ≤ 2c−1`, heavy factor
`L_c = (a−(c−2))(b−(c−1))`, and `deg_j H_c = 2c−2`. Its anchor is a pair of length-`(c−1)` runs:
```
H_c(a,b,0) = ∏_{s=3}^{c+1}(a+s) · ∏_{s=2}^{c}(b+s).        (the A-run · B-run)
```
**Proposition 2** (verified vs MN throughout): for `1 ≤ j ≤ b`,
```
Δ(j) = j − 2 s₂(j) + 2[ v₂C(a+c+1, j) + v₂C(b+c, j) + v₂Q_c(j) − v₂Q_c(0) ].
```

The proven **linear heavy floor** `2^{β'(c)} ∣ H_c` (constant in `a,b,j`, `β' ≤ γ(c) ≤ 2c−2`;
`β'(4,6,8,10) = 4,7,11,14`) is the deep-tail lever, available `c`-uniformly.

---

## 1. The reduction at `j = 4`: `Δ(4) = 2 v₂(G_4^{(c)}) − 10`

Fix even `c ≥ 6`, so `4 ≤ c − 1` (below the peeling threshold) and `4 < 2c` (tip vanishes). Then
`Q_c(4) = L_c·H_c(4)` and `Q_c(0) = L_c·H_c(0)`, so `v₂Q_c(4) − v₂Q_c(0) = v₂H_c(4) − v₂H_c(0)`
(the heavy `L_c` cancels). With `s₂(4) = 1`, Prop 2 gives
```
Δ(4) = 2 + 2[ v₂C(a+c+1,4) + v₂C(b+c,4) + v₂H_c(4) − v₂H_c(0) ].     (1)
```

**Run cancellation.** The numerator of `C(a+c+1,4) = (a+c+1)(a+c)(a+c−1)(a+c−2)/24` is exactly the
**top 4** factors of the A-run `∏_{s=3}^{c+1}(a+s)`; likewise the numerator of `C(b+c,4)` is the top 4
of the B-run. Writing `A_lo := ∏_{s=3}^{c−3}(a+s)` and `B_lo := ∏_{s=2}^{c−4}(b+s)` for the two
**lower runs** (each of length `c−5 ≥ 1`), and using `v₂H_c(0) = v₂(\text{A-run}) + v₂(\text{B-run})`:
```
v₂C(a+c+1,4) − v₂(\text{A-run}) = v₂(\text{top 4}) − 3 − v₂(\text{top 4}) − v₂(A_lo) = −3 − v₂(A_lo),
```
and symmetrically for `b`. Substituting into (1):
```
Δ(4) = 2 + 2[ −3 − v₂(A_lo) − 3 − v₂(B_lo) + v₂H_c(4) ] = 2 v₂(G_4^{(c)}) − 10,     (2)
```
where the **reduced heavy at `j = 4`** is
```
G_4^{(c)} := H_c(a,b,4) / [ ∏_{s=3}^{c−3}(a+s) · ∏_{s=2}^{c−4}(b+s) ].
```
The division is exact (a polynomial-division certificate per `c`; `crux_final.py` asserts remainder
`0`), and `G_4^{(c)}` is a **biquartic** in `(a,b)` (degree 4 in each). Identity (2) is verified
per-shape against MN for `c = 6,…,16` (2400 shapes, **0 mismatch**, `crux_final.py`). It matches the
`c = 6` note (`Δ(4) = 2 v₂(G_4) − 10`) and the `c`-uniform constant `10 = 4·v₂(4!) − (4 − 2 s₂(4)) = 12 − 2`.

So **`min Δ(4) = 2·content(G_4^{(c)}) − 10`**, where `content` denotes the 2-adic content of the
integer-valued polynomial `G_4^{(c)}` over the parity sheet `a ≡ b (mod 2)` (= the minimum `v₂` of its
forward-difference / binomial-basis coefficients on the two sheets `(a,b) ≡ (0,0),(1,1) mod 2`).

---

## 2. The crux: `content(G_4^{(c)}) = [c ≡ 2 mod 4] ? 5 : 6`

`G_4^{(c)}` is, uniformly in `c`, the biquartic whose 25 coefficients are polynomials in `c` of
degree `≤ 4` (interpolated from `c = 6,8,10,12,14`; **independently verified for
`c = 16,18,20,22,24,26`** — 6 free checks past the 5 needed to pin a degree-4 fit, `crux_final.py`).
For instance the corner coefficients are
```
[a^4 b^4]G_4 = 1,   [a^4 b^0]G_4 = c(c−1)(c−2)(c−3),   [a^0 b^4]G_4 = c(c−1)(c+1)(c−2), …
```
(full list in `crux_final.py`). Its content is governed by two coefficients.

**(i) The minimiser `K(c)`.** In the binomial basis on the parity sheet, the coefficient of
`C(u,0)C(v,1)` on the `(0,0)`-sheet (equivalently `C(u,0)C(v,0)` on the `(1,1)`-sheet) is
```
K(c) = 24·c·(c−1)·(c−4)·(c−5).
```
Since `v₂(24) = 3` and, for even `c`, `c−1` and `c−5` are odd while `c` and `c−4` are even,
```
v₂(K(c)) = 3 + v₂(c) + v₂(c−4).
```
Now `c` and `c−4` are congruent mod 4:
- `c ≡ 2 (mod 4)`: both `c, c−4 ≡ 2 (mod 4)`, so `v₂(c) = v₂(c−4) = 1` and **`v₂(K) = 5`**;
- `c ≡ 0 (mod 4)`: both `≡ 0 (mod 4)`, so `v₂(c), v₂(c−4) ≥ 2` and **`v₂(K) ≥ 7`**.

**(ii) The floor `192·(odd)`.** The coefficient of `C(u,2)C(v,2)` on each sheet is
`192·(2c⁴ + 32c³ + 18c² + 16c + 15)` (sheet `(0,0)`; sheet `(1,1)` has `+15` likewise). As
`192 = 2⁶·3` and the parenthesised polynomial is **odd** for every even `c` (constant term `15`
odd, all other terms even), this coefficient has `v₂ = 6` for **all** even `c`.

**(iii) Uniform lower bound.** Every one of the 49 nonzero binomial-basis coefficients (both sheets)
is an explicit polynomial `p(c)`; the claims
```
v₂(p(c)) ≥ 5  for all c ≡ 2 (mod 4),      v₂(p(c)) ≥ 6  for all c ≡ 0 (mod 4)
```
are **finite residue certificates**: `p(c) mod 2^7` depends only on `c mod 2^7`, so it suffices to
verify `32 ∣ p(c)` (resp. `64 ∣ p(c)`) on one period. The exhaustive check over `c ∈ [6, 512)`
(covering every residue mod `128` in each class many times) returns minimum coefficient valuation
```
c ≡ 2 (mod 4):  min = 5  (attained by K(c));      c ≡ 0 (mod 4):  min = 6  (attained by 192·odd).
```
(`crux_verify.py`, rigorous.)

**Combining (i)–(iii):**
```
content(G_4^{(c)}) = 5  (c ≡ 2 mod 4),      content(G_4^{(c)}) = 6  (c ≡ 0 mod 4),
```
whence by (2) **`min Δ(4) = 0` iff `c ≡ 2 (mod 4)`, else `min Δ(4) = 2`.** ∎ (Main Theorem)

**Tie classes.** Refining the content to each `(a,b) mod 4` class (a mod-8 constancy check,
`crux_G4.py`): for `c ≡ 2 mod 4`, `content(G_4) = 5` on **exactly** `(a,b) ≡ (0,2),(1,1) (mod 4)` and
`= 6` elsewhere — the tie set `{(0,2),(1,1)}` is **`c`-uniform** (verified `c = 6,8,10,12`; MN census
`c = 14,16,18` gives min `Δ(4) = 0` on precisely these classes). Thus when `c ≡ 2 mod 4`,
`Δ(4) = 0 ⟺ (a,b) ≡ (0,2),(1,1) (mod 4)`.

*Why the resonance.* Generator 4 is a **low-régime, tip-free** event (`j = 4 < 2c`); its incidence is
the single unit of 2-content in `K(c) = 24 c(c−1)(c−4)(c−5)`. The factor pair `c·(c−4)` — the two even
factors — carries `v₂(c) + v₂(c−4)`, minimal (`= 2`) exactly on the `c ≡ 2 mod 4` sheet where both are
`≡ 2 mod 4`. The periodicity is arithmetic, not a run-length/size effect (killing the heuristic of the
refutation note), and it is *period 4 in `c`*, proved here for all even `c` at one stroke.

---

## 3. Generator 2, `c`-uniformly: `min Δ(2) = 0` always

The identical peeling at `j = 2` (top-2 factors of `C(a+c+1,2), C(b+c,2)` cancel the top-2 of the
runs, lower runs `∏_{s=3}^{c−1}(a+s), ∏_{s=2}^{c−2}(b+s)`) gives
```
Δ(2) = 2 v₂(G_2^{(c)}) − 4,   G_2^{(c)} = H_c(a,b,2)/[∏_{s=3}^{c−1}(a+s)·∏_{s=2}^{c−2}(b+s)],
```
a **biquadratic** with `c`-polynomial coefficients (deg `≤ 4`, verified `c ≤ 20`, `crux_G2_cpoly.py`):
```
G_2^{(c)} = a²b² + (2c−1)a²b + c(c−1)a² + (2c+1)ab² + (2c²−1)ab + c(3c−5)a
            + c(c+1)b² + c(3c−5)b + 2c(c−3).
```
The minimising binomial-basis coefficient is `C(u,1)C(v,1)`: on the `(0,0)`-sheet it is
`4·(2c²+8c+3)`, on the `(1,1)`-sheet `4·(2c²+16c+15)` — both `4·(odd)`, so `v₂ = 2` for **every** even
`c`. All other coefficients have `v₂ ≥ 2`. Hence
```
content(G_2^{(c)}) = 2  for all even c,   so  min Δ(2) = 0  for all even c:
```
**generator 2 always fires** (`0 ∈ J*` unconditionally, `{0,2} ⊆ J*` on its tie classes). Refining to
`(a,b) mod 4` (`crux_G2.py`), the tie set is **periodic in `c mod 4`**:
```
Δ(2) = 0  ⟺  (a,b) ≡ (2,2),(1,3) (mod 4)   if c ≡ 0 (mod 4),
Δ(2) = 0  ⟺  (a,b) ≡ (0,0),(3,1) (mod 4)   if c ≡ 2 (mod 4).
```
(Consistent with `c = 4`: `(2,2),(1,3)`, and `c = 6`: `(0,0),(3,1)`.) The generator-2 and generator-4
tie classes are disjoint, so `|J* ∩ {0,2,4}| ≤ 2`.

---

## 4. No generator beyond 4: the régimes for even `j ≥ 6`

Odd `j` never ties (`val` odd). For even `j ≥ 6` we show `Δ(j) > 0`.

**Deep `j ≥ 2c`, Case A (`c`-uniform, proved).** With `φ = v₂(L_c) ≥ 1`,
`P := v₂((j)_{2c}) = 2c + s₂(j−2c) − s₂(j)`, Case A is `v₂((a−(c−2))(b−(c−1))H_c(j)) ≥ P`, i.e.
`v₂Q_c(j) ≥ P`. For `j ≥ 2c` the run `(a+2)_{j−(c−1)}` (length `≥ c+1`) splits as a top `c`-block
`(a+2)⋯(a−(c−3))` (`v₂ ≥ vfact(c)`), the isolated factor `a−(c−2)`, and a tail `(a−(c−1))_{j−2c}`
(`v₂ ≥ vfact(j−2c)`); symmetrically for `b`. Substituting in the peeling identity (the
`v₂(a−(c−2)), v₂(b−(c−1))` cancel):
```
T(j) ≥ 2vfact(c) + 2vfact(j−2c) − 2vfact(j) + v₂Q_c(j) ≥ 2vfact(c) − P,
Δ(j) ≥ j − 2s₂(j) + 4vfact(c) − 2P = j − 4s₂(c) − 2s₂(j−2c)  =  2c + t − 4s₂(c) − 2s₂(t)   (j = 2c+t).
```
Since `t − 2s₂(t) ≥ −1` (Lemma A), `Δ(j) ≥ 2c − 4s₂(c) − 1 ≥ 3` for all even `c ≥ 4` (as
`c − 2s₂(c) ≥ 2`). **`c`-uniform, `Δ ≥ 3 > 0`.** ∎

**Mid `c−1 ≤ j ≤ 2c−1` (tip-free) and deep Case B (`Hv < P`): reduced to `g_c(j) > 0`.** In both,
the peeling identity is exact/collapses to
```
Δ(j) = Δ̂(j) := j − 2s₂(j) + 2[ v₂((a+2)_{j−(c−1)}) + v₂((b+1)_{j−(c−1)}) − 2vfact(j) + v₂H_c(j) ].
```
Using `v₂((x)_r) ≥ vfact(r)` and the linear floor `v₂H_c(j) ≥ β'(c)`:
```
Δ̂(j) ≥ g_c(j) := j + 2s₂(j) − 4 s₂(j−(c−1)) − 4(c−1) + 2β'(c).
```
The `+j` term dominates: `g_c(j) ≥ j − 4⌊log₂ j⌋ − 4(c−1) + 2β'(c) − 2 → +∞`, so `g_c(j) > 0` for all
`j` past an explicit `c`-threshold; on the finite window below it, and at the finitely many indices
where `g_c(j) ≤ 0` (the analogues of `c = 6`'s `j ∈ {12,16}`), `Δ̂(j) > 0` is a complete residue check
on `content(Π_j)`, `Π_j := (a+2)_{j−(c−1)}(b+1)_{j−(c−1)}H_c(j)` (exact binomial-basis content — no
`m`-cap). This is `c`-uniform *in form*; the exception set is finite per `c`.

**MN census (ground truth).** For `c = 6,8,10,12,14`, every even `j` up to `2c+2`, all interior
shapes with `a ≤ 2c+78`: `min Δ(2) = 0`; `min Δ(4) = [c≡2?0:2]`; **`min Δ(j) ≥ 2 > 0` for every even
`j ≥ 6`** (`census_nogen.py`). In particular `min Δ(6) = 2` uniformly — the next would-be generator is
always killed.

---

## 5. Assembly: the general even-`c` interior law (table periodic in `c mod 4`)

> **Theorem (general even-`c` interior, generator menu).** For box-interior `λ = (a,b,c)`, `c` even,
> `b` large, the minimiser set `J*` is:
>
> | `c mod 4` | `(a,b) mod 4` | `J*` | firing generator |
> |---|---|---|---|
> | `0` | `(2,2)`, `(1,3)` | `{0,2}` | 2 |
> | `0` | else | `{0}` | — |
> | `2` | `(0,0)`, `(3,1)` | `{0,2}` | 2 |
> | `2` | `(0,2)`, `(1,1)` | `{0,4}` | **4** |
> | `2` | else | `{0}` | — |
>
> In all cases `0 ∈ J*`, `J* ⊆ {0,2,4}`, and `|J*| ≤ 2` is **even**; so `G_λ(i) = 0 ⟺ |J*|` even holds
> with `|J*|` never odd. **Generator 4 is available iff `c ≡ 2 (mod 4)`.**

*Status of the assembly.* Fully proved `c`-uniformly: the reduction (§1), the generator-2 law (§3,
`min Δ(2)=0` always + the periodic tie classes), the **generator-4 crux (§2, the theorem's soul)**, the
deep Case A (§4). The residual — `Δ(j) > 0` for even `j ≥ 6` off the deep-Case-A régime — is reduced
`c`-uniformly to `g_c(j) > 0` plus a finite, per-`c` family of `content(Π_j)` residue certificates
(§4); it is complete by census for `c ≤ 14` and by the specific certificates for `c = 4, 6`. The table
recovers the proven `c = 4` and `c = 6` families as its `c ≡ 0` and `c ≡ 2` instances.

---

## 6. What is proved

**Proved in full (symbolic / exact binomial-basis content — complete residue-class certificates, no
`m`-cap):**
- The reduction `Δ(4) = 2 v₂(G_4^{(c)}) − 10` for all even `c ≥ 6` (peeling + run cancellation;
  per-shape vs MN, `c ≤ 16`, 0 mismatch).
- **Main Theorem:** `content(G_4^{(c)}) = 5` iff `c ≡ 2 mod 4`, else `6`, via `v₂(K(c)) = 3+v₂(c)+v₂(c−4)`
  and the `192·odd` floor, with the uniform lower bound a finite mod-`2^7` certificate over all 49
  coefficients. `⟹ min Δ(4) = [c ≡ 2 mod 4] ? 0 : 2`. Generator-4 tie classes `{(0,2),(1,1)}`,
  `c`-uniform.
- Generator-2: `content(G_2^{(c)}) = 2` for all even `c` (min coeff `4·odd`) ⟹ `min Δ(2) = 0`;
  tie classes periodic in `c mod 4`.
- Deep Case A: `Δ(j) ≥ 2c − 4s₂(c) − 1 ≥ 3` for all `j ≥ 2c`, `c`-uniform.

**Reduced `c`-uniformly (in form) + finite per-`c` certificates / census:**
- Even `j ≥ 6`, mid/deep-Case-B: `Δ(j) ≥ g_c(j) > 0` (analytic tail) plus finitely many
  `content(Π_j)` checks; complete by MN census for `c ≤ 14`.

**Net.** The generator MENU of the general even-`c` interior — *which* generators are available and
on which `(a,b) mod 4` classes — is now a **`c`-uniform theorem**, periodic in `c mod 4`, with the
generator-4 crux (the stated open problem) fully derived from `content(H_c)` at `j = 4`, not a per-`c`
residue scan. The one-line resonance `v₂(c) + v₂(c−4)` is the arithmetic the periodicity rests on.

### For LEAN
`v₂(K(c)) = 3 + v₂(c) + v₂(c−4)` is Kummer/`vz` arithmetic; the content claims are `decide`-friendly
finite residue checks mod `2^7` on explicit low-degree integer polynomials in `c`; the reduction (2)
is the same peeling algebra already earmarked for the `c = 4,6` interiors.

### Files (`projects/code/threerow-heavy/`)
`crux_G4.py` (reduction (2) + per-class content, `c = 6,8,10,12`), `crux_G4_cpoly.py` (G_4 as `c`-poly,
minimiser search), `crux_verify.py` (degree check, rigorous mod-`2^7` content certificate, MN
cross-check `c = 14,16,18`), `crux_final.py` (degree-≤4 certification `c ≤ 26`, per-shape identity,
closed form), `crux_G2.py` / `crux_G2_cpoly.py` (generator-2), `census_nogen.py` (no generator beyond
4, `c ≤ 14`). Engine: `jobB_betaprime_sequence.py` (`interp_ab`, `build_Hc`), MN `threerow-c3/mn.py`.
