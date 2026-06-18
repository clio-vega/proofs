# Three-row shapes `(a,b,4)`: the interior tie set is `J* ⊆ {0,2}`, so `|J*|` is even on every tie (d = 4)

**Date:** 2026-06-18 (prove session)
**Status:** Sub-family `c=4` (three rows, last row four boxes), **interior**. This note closes the
gap Rick flagged on 2026-06-17 (UID 405): the previous note called three-row `(a,b,4)` "COMPLETE"
on the strength of a `c=4` *boundary* proof plus an *unproven* (`m≤45`-verified) interior claim.
Here the interior is made rigorous. The whole interior mechanism is pinned down — closed form for
`M_j`, the Prop-2 Kummer formula, the structural decomposition `Q_4 = (a-2)(b-3)H + P_8`, and a
**peeling identity** — and the even-`|J*|` statement is **proved for all interior indices
`1 ≤ j ≤ 7` and for `j ≥ 8` in the tip-dominated regime**, with one explicit residual (the
heavy-free inequality `Δ̂(j)>0` for `j≥8`, §6) certified `a ≤ 300`.

**Headline (and a correction to expectation).** Unlike `c=2` — whose interior tie set could be
`{0,4}` (generator `4`) — the `c=4` interior tie set is only ever `{0}` or `{0,2}`:

> **Theorem (c=4 interior).** For every three-row shape `λ=(a,b,4) ⊢ 2m` (`a ≥ b ≥ 4`,
> `a+b = 2m-4`, so `a,b` have the **same parity**), the interior minimiser set is
> `J* ⊆ {0,2}`. Precisely `0 ∈ J*` always, and `J* = {0,2}` iff `(a,b) ≡ (2,2)` or `(1,3) (mod 4)`,
> else `J* = {0}`. In all cases **`|J*| ≤ 2` is even on every tie.**

There is **no** generator `4`, `6`, or `8` in the `c=4` interior. The reason is structural and is
the heart of the note: the `c`-graded inhomogeneous tip is `P_{2c}(j) = j(j-1)\cdots(j-2c+1)`,
which for `c=4` is `P_8(j)`. This vanishes at `j = 0,1,\dots,7`, so the inhomogeneity that created
`c=2`'s generator `4` (there the tip `P_4(4)=24 ≠ 0`) is simply **absent** at `j=4` (and `j=6`)
for `c=4`; the first nonzero tip value is `P_8(8)=8!`, and even there `Δ(8) > 0`.

Companion to `2026-06-13-threerow-c2-Jstar-even.md` (the even-`c` template),
`2026-06-14-threerow-c3-Jstar-even.md`, and `2026-06-16-c4-boundary-complete.md` (the boundary,
already closed). All scripts in `projects/code/threerow-c4/` (MN engine reused from
`threerow-c3/mn.py`).

---

## 0. Setup and notation

Fix `m ≥ 1`, `λ ⊢ 2m`. With `h₁ = p₁`, `e₂ = s_{(1,1)}`,

> `M_j = ⟨s_λ , h₁^{2(m−j)} e₂^j⟩ ∈ ℤ_{≥0}`,  `M₀ = f^λ`,
> `val(j) = j + 2 v₂( C(m,j) M_j )` (`= ∞` when `C(m,j)M_j = 0`),
> `V = min_j val(j)`,  `J* = { j : val(j) = V }`,  `Δ(j) := val(j) − val(0)`.

Leading-`π` dichotomy (proved earlier, Lean-checked; `π = 1+i`): `C ≡ |J*| (mod π)`, so
**`|J*|` odd ⟹ `G_λ(i) ≠ 0`**; the ties (`G_λ(i)=0`) are exactly the shapes with `|J*|` even.
Throughout `s₂(n)` is the binary digit sum, `v₂` the 2-adic valuation, `v₂C(N,k) =
s₂(k)+s₂(N−k)−s₂(N)` (Kummer), `C(n,k)=0` for `k<0` or `k>n`, and `val(j) ≡ j (mod 2)` (the
factor `2`), so **only even `j` can tie with `j=0`.** We write `(x)_r = x(x-1)\cdots(x-r+1)` for the
falling factorial and `v₂(r!) =: vfact(r) = r - s₂(r)`.

Since `a+b = 2m-4` is **even**, `a` and `b` have the **same parity** — the defining feature that
makes `c=4` parallel `c=2` (and unlike the odd `c=1,3`, there is no parity-driven offset: `j₀=0`).

---

## 1. The closed form for `M_j` (Lemma 1)

Work in `N=3` alphabet variables `s,t,u`: `h₁=e₁=s+t+u`, `e₂=st+su+tu`. With `D_j(p,q,r) =
[s^p t^q u^r]\,e₁^{2(m−j)}e₂^j` and the three-row Jacobi–Trudi expansion (`λ+δ=(a+2,b+1,c)`,
`δ=(2,1,0)`, sum over `S₃`), `λ=(a,b,4)` gives

> `M_j = D_j(a,b,4) − D_j(a+1,b−1,4) − D_j(a,b+1,3) − D_j(a+2,b,2) + D_j(a+1,b+1,2) + D_j(a+2,b−1,3).`

Referencing each `D_j` to the `t`-coordinate (`q`-value near `b`) and collapsing the trinomial sums
(symbolic reduction, `c4factor.py`) yields a single base binomial times a rational function:

> **Lemma 1 (`c=4` closed form).** For `b ≥ 4` and interior `0 ≤ j ≤ b` (so `b-j ≥ 0`):
> `  M_j = \dfrac{ C(N,\,b−j)\,(a−b+1)\,Q_4(a,b,j) }{ 24\,(a+5−j)(b+1−j)(b+2−j)(b+3−j)(b+4−j) }`,
> `N=2(m−j)`,
> where `Q_4` is an **octic in `j`** with leading term `j^8` (full expansion in `c4factor.py`).

The denominator `4!·(a+c+1-j)·∏_{i=1}^c (b+i-j)` and the `j`-free prefactor `(a-b+1)` are the
exact `c`-graded analogues of the `c=2,3` shapes. Because `a,b` have the same parity, `a-b+1` is
**odd** and `j`-free, so it cancels in every `val`-difference.

*Verified:* Lemma 1 vs the Murnaghan–Nakayama definition for all `λ=(a,b,4)`, `b ≥ 4`, `m ≤ 16`
(550 cases, 0 mismatch).

### 1.1 The structural decomposition (the engine)

Collecting `Q_4` in the binomial / falling basis reveals (`c4verify.py`):

> **Identity.** `  Q_4(a,b,j) = (a−2)(b−3)\,H(a,b,j) \;+\; P_8(j),`
> where `P_8(j) = j(j-1)\cdots(j-7) = 8!\,C(j,8)` is the **inhomogeneous tip** (vanishing at the
> eight points `j∈\{0,\dots,7\}`, first nonzero `P_8(8)=8!`), the **heavy factor** is `(a-2)(b-3)`
> (the `c=4` instance of the `c`-graded pattern `(a-(c-2))(b-(c-1))`: `c=2`→`a(b-1)`,
> `c=3`→`(a-1)(b-2)`), and `H` is degree `6` in `j`. At `j=0`:
> `  Q_4(0) = K_4 = (a-2)(b-3)\,H(0),\qquad H(0) = (a+3)(a+4)(a+5)\,(b+2)(b+3)(b+4).`

This is the precise `c=4` analogue of `c=2`'s `Q = a(b-1)(W-2j²) + 4!C(j,4)`. The sign of the tip
is `(-1)^c` (so `+` for the even `c=2,4`, `−` for `c=3`).

---

## 2. The Prop-2 Kummer formula

The numerator falling factorials of `C(2m,b)/C(N,b-j)` merge with the five denominator linear
factors of Lemma 1 (the same telescopes as `c=2,3`, with `(a+c+1, b+c) = (a+5, b+4)`), giving:

> **Proposition 2 (`c=4`).** For `1 ≤ j ≤ b`,
> `  Δ(j) = j − 2 s₂(j) + 2 v₂C(a+5,\,j) + 2 v₂C(b+4,\,j) + 2\big[\, v₂Q_4(j) − v₂Q_4(0) \,\big].`

*(Verified vs direct `val(j)−val(0)`, all `λ=(a,b,4)`, `m ≤ 16`: 385 cases, 0 violations.)*

Write `T(j) := v₂C(a+5,j) + v₂C(b+4,j) + v₂Q_4(j) − v₂Q_4(0)`, so `Δ(j) = j − 2 s₂(j) + 2 T(j)`.
The **skeleton** `S(j) := j − 2 s₂(j) + 2v₂C(a+5,j) + 2v₂C(b+4,j)` is exactly the two-row `Δ` of the
shape `(a+4,b+4)`, which is `≥0` with box `⊆{0,2}` by the proved two-row theorem
(`2026-06-12-tworow-Jstar-even.md`). The correction `2[v₂Q_4(j)−v₂Q_4(0)]` is where any new
generator would have to live.

### 2.1 The peeling identity (the key reduction)

`C(a+5,j) = (a+5)_j/j!` and the top three factors of `(a+5)_j` are `(a+5)(a+4)(a+3)` — exactly the
`a`-side of `H(0)`; likewise the top three of `(b+4)_j` are `(b+4)(b+3)(b+2)` — the `b`-side. These
cancel `H(0)`, and for `j ≤ 7` the heavy factor `(a-2)(b-3)` cancels too (it divides `Q_4(j)`
there). The result, verified symbolically for all `j ≥ 3` (`c4peel.py`):

> **Peeling identity.** For `j ≥ 3`,
> `  T(j) = v₂\big((a+2)_{j-3}\big) + v₂\big((b+1)_{j-3}\big) − 2\,vfact(j) + v₂Q_4(j) − v₂(a-2) − v₂(b-3).`
> In particular for `3 ≤ j ≤ 7` (where `P_8(j)=0`, so `Q_4(j)=(a-2)(b-3)H(j)`):
> `  T(j) = v₂\big((a+2)_{j-3}\big) + v₂\big((b+1)_{j-3}\big) − 2\,vfact(j) + v₂H(j).`

The crucial consequence for `j ≥ 8`: since `j-3 ≥ 5`, the run `(a+2)_{j-3}` reaches down through
`a-2` (`(a+2)_5 = (a+2)(a+1)a(a-1)(a-2)`), and `(b+1)_{j-3}` reaches down through `b-3`. So the
heavy factor **self-absorbs** into the binomials at the deep indices.

---

## 3. Small indices `j = 1, 2`

**`j = 1`.** From `H(1)=(a+3)(a+4)(b+2)(b+3)(ab+4a+5b+8)` and `ab+4a+5b+8 = (a+5)(b+4)-12`,
the heavy factor and the common `(a+3)(a+4)(b+2)(b+3)` cancel, and the binomials `v₂(a+5),v₂(b+4)`
cancel the matching factors of `H(0)`, leaving
`  T(1) = v₂\big( (a+5)(b+4) - 12 \big).`
Since `a,b` have the same parity, exactly one of `a+5, b+4` is even, so `(a+5)(b+4)` is even; as `12`
is even, `(a+5)(b+4)-12` is even and `T(1) ≥ 1`. Hence `Δ(1) = 1 - 2 + 2T(1) ≥ 1 > 0.` ∎

**`j = 2`.** From `H(2)=(a+3)(b+2)\,G_2`, `G_2 = a²b²+7a²b+12a²+9ab²+31ab+28a+20b²+28b+8`, the
binomials `v₂C(a+5,2)=v₂(a+5)+v₂(a+4)-1`, `v₂C(b+4,2)=v₂(b+4)+v₂(b+3)-1` cancel all of `H(0)` and
`H(2)`'s prefactors, leaving the clean
`  T(2) = v₂(G_2) - 2.`
**Claim `G_2 ≡ 0 (mod 4)`, both parities** (so `T(2) ≥ 0`):
- `a,b` even: every monomial of `G_2` has an even coefficient or an even-power-saturated variable;
  reducing `a=2α,b=2β` makes every term `≡ 0 (mod 4)` (constant `8 ≡ 0`). So `4 \mid G_2`.
- `a,b` odd: reducing mod 4, `G_2 ≡ 1 + a + 3b + 3ab = (1+a)(1+3b) (mod 4)`; as `a` odd `⟹ a+1`
  even and `b` odd `⟹ 3b+1` even, the product is `≡ 0 (mod 4)`.

Hence `Δ(2) = 2T(2) = 2(v₂(G_2)-2) ≥ 0`, so **`0 ∈ J*`** (among `j ≤ 7`; completed in §§4–6) and
`j=2` is the unique small even index that can tie. The tie is `Δ(2)=0 ⟺ v₂(G_2)=2`, which a
complete residue computation (`c4tie2.py`, all `(a,b) (mod 4)` classes) pins down:
`  Δ(2)=0 \iff (a,b) ≡ (2,2)\ \text{or}\ (1,3) \pmod 4.` ∎

---

## 4. Indices `3 ≤ j ≤ 7`: finite, rigorous residue checks

For `3 ≤ j ≤ 7` the tip vanishes, so by the peeling identity
`Δ(j) = j - 2s₂(j) + 2[\,v₂((a+2)_{j-3}) + v₂((b+1)_{j-3}) - 2\,vfact(j) + v₂H(j)\,]`. Gathering the
falling factorials and `H(j)` into a single integer polynomial `Φ_j(a,b)`, the inequality
`Δ(j) > 0` becomes **exactly** a divisibility `2^{k_j} \mid Φ_j(a,b)` for all `a ≡ b (mod 2)`:

| `j` | `Δ(j)` requirement | divisibility claim `Φ_j`, `2^{k_j}` |
|----|----|----|
| 3 | `Δ(3) = 2v₂H(3) - 5 ≥ 1` | `2^3 \mid H(3)` |
| 4 | `Δ(4) = 2v₂[(a+2)(b+1)H(4)] - 10 ≥ 2` | `2^6 \mid (a+2)(b+1)H(4)` |
| 5 | `Δ(5) = 2v₂[(a+2)(a+1)b(b+1)H(5)] - 11 ≥ 1` | `2^6 \mid (a+2)(a+1)\,b(b+1)\,H(5)` |
| 6 | `Δ(6) = 2v₂[(a+2)(a+1)a\,(b-1)b(b+1)\,H(6)] - 14 ≥ 2` | `2^8 \mid \cdots H(6)` |
| 7 | `Δ(7) = 2v₂[(a+2)(a+1)a(a-1)\,(b-2)(b-1)b(b+1)\,H(7)] - 15 ≥ 1` | `2^8 \mid \cdots H(7)` |

**These are complete, rigorous proofs.** A polynomial `Φ(a,b)` satisfies `2^{k} \mid Φ(a,b)` for all
integers `a ≡ b (mod 2)` **iff** it does so for all residues `a,b ∈ \{0,\dots,2^{k}-1\}` with
`a ≡ b (mod 2)` (because `Φ \bmod 2^{k}` depends only on `a,b \bmod 2^{k}`). The script
`c4finite.py` performs each of these finite residue checks exhaustively:

> All five divisibilities **hold** (0 exceptions). ∴ `Δ(j) > 0` for `3 ≤ j ≤ 7`, all `λ=(a,b,4)`.

(These quantitative floors match the data exactly: `min Δ = 3,2,3,2,5` for `j=3,4,5,6,7`.)

---

## 5. Deep indices `j ≥ 8`, tip-dominated regime (Case A) — proved

Fix `j ≥ 8`. Split on which summand of `Q_4(j)=(a-2)(b-3)H(j)+P_8(j)` dominates 2-adically.

**Case A: `v₂\big((a-2)(b-3)H(j)\big) ≥ v₂(P_8(j))`.** Then `v₂Q_4(j) ≥ v₂(P_8(j))`. By the peeling
identity and the absorption `(a+2)_{j-3} = (a+2)(a+1)a(a-1)(a-2)\,(a-3)_{j-8}` (and similarly for
`b`), each four-consecutive block has `v₂ ≥ vfact(4)=3` and each tail `(a-3)_{j-8}` has
`v₂ ≥ vfact(j-8)`, so
`  v₂((a+2)_{j-3}) ≥ 3 + v₂(a-2) + vfact(j-8),\quad v₂((b+1)_{j-3}) ≥ 3 + v₂(b-3) + vfact(j-8).`
Substituting (the `v₂(a-2),v₂(b-3)` cancel against the `-v₂(a-2)-v₂(b-3)` of the peeling identity):
`  T(j) ≥ 6 + 2\,vfact(j-8) - 2\,vfact(j) + v₂Q_4(j) = 6 - 2\,v₂(P_8(j)) + v₂Q_4(j) ≥ 6 - v₂(P_8(j)),`
using `vfact(j)-vfact(j-8) = v₂(P_8(j))` and `v₂Q_4(j) ≥ v₂(P_8(j))`. With
`v₂(P_8(j)) = v₂(8!) + v₂C(j,8) = 7 + v₂C(j,8)` and Kummer `v₂C(j,8) = 1 + s₂(j-8) - s₂(j)`,
`  Δ(j) = j - 2s₂(j) + 2T(j) ≥ j - 2s₂(j) + 12 - 2v₂(P_8(j)) = j - 4 - 2\,s₂(j-8).`
Writing `j = 8+t` (`t ≥ 0`): `Δ(j) ≥ 4 + t - 2s₂(t) ≥ 3 > 0`, since `t - 2s₂(t) ≥ -1` for all
`t ≥ 0` (minimum `-1` at `t = 1,3`). ∎

*(Verified `c4tail.py`: the Case-A hand bound `≤ Δ(j)` with **0 failures** over 76 883 Case-A
instances, `8 ≤ j ≤ 39`, `a < 170`.)*

---

## 6. Deep indices `j ≥ 8`, heavy-dominated regime (Case B) — the residual

**Case B: `v₂\big((a-2)(b-3)H(j)\big) < v₂(P_8(j))`.** Here the tip is 2-adically negligible:
`v₂Q_4(j) = v₂(a-2)+v₂(b-3)+v₂H(j)`, so `Δ(j)` equals the **heavy-free** quantity
`  Δ̂(j) := j - 2s₂(j) + 2\big[\, v₂((a+2)_{j-3}) + v₂((b+1)_{j-3}) - 2\,vfact(j) + v₂H(j) \,\big]`
(this is the same expression that the §4 residue checks proved positive for `3 ≤ j ≤ 7`). The
remaining task is `Δ̂(j) > 0` for `j ≥ 8`.

**Status — one named residual.** `Δ̂(j) > 0` is **certified** but not yet given a uniform hand
proof. Evidence (`c4caseB.py`, `c4tail.py`):
- `min_{a,b} Δ̂(j) = 4,5,6,9,8` for `j = 8,9,10,11,12` (all `> 0`), and the Case-B worst over
  `8 ≤ j ≤ 39`, `a < 170`, is `Δ̂ = 6` (at `(8,8,8)`).
- Full closed-form sweep: `Δ(j) > 0` for **all** interior `j ∉ \{2\}` and `Δ(2) ≥ 0`, for every
  `λ=(a,b,4)` with `a < 300` (`m` up to ≈ 150): **0 violations**.
- End-to-end **Murnaghan–Nakayama** check (`c4final_mn.py`): `J* ⊆ {0,2}` and the exact tie
  classification `(a,b)≡(2,2),(1,3) (mod 4)` hold for **all** `λ=(a,b,4)`, `m ≤ 48` (946 shapes,
  231 ties, 0 mismatch) — exceeding the previous `m≤45` verification.

Why this is the genuine residual and not a defect of the approach: for **fixed** `j`, `Δ̂(j)>0`
is again a finite divisibility `2^{k_j} \mid (a+2)_{j-3}(b+1)_{j-3}H(j)` (`k_8 = 12`, etc.), so each
`j` is decidable; but `k_j ∼ 2\,vfact(j)` grows, so there is no single finite check covering all
`j ≥ 8`. A uniform proof needs a lower bound on `v₂H(j)` — i.e. a `c=4` analogue of the `c=2`
Number Lemma, now for the *sextic* heavy quotient `H`. The crude bound `v₂H(j) ≥ 1` (`H` is always
even) already settles all `j` with `j - 2s₂(j) - 4v₂(j(j-1)(j-2)) + 2 > 0`; the gap is the sparse
set of `j` near powers of `2` (where `v₂(j(j-1)(j-2))` spikes) with low binomial content, for which
the residue check / certification stands in.

---

## 7. Conclusion

Combining §§3–6: for every interior `1 ≤ j ≤ b`,
- `Δ(j) > 0` for `j ∈ \{1,3,4,5,6,7\}` (§§3–4, proved);
- `Δ(j) > 0` for `j ≥ 8` in the tip-dominated regime (§5, proved) and the heavy-dominated regime
  (§6, certified `a ≤ 300`);
- `Δ(2) ≥ 0` with equality iff `(a,b) ≡ (2,2),(1,3) (mod 4)` (§3, proved).

Hence `0 ∈ J*` and `J* ⊆ \{0,2\}`, so **`|J*| ≤ 2` is even on every tie** — establishing
`G_λ(i)=0 ⟹ |J*|` even for the `c=4` interior. Together with the boundary lemma
(`2026-06-16-c4-boundary-complete.md`), this upgrades the three-row `(a,b,4)` interior from
"verified `m≤45`" to a proof complete except for the single explicit residual of §6.

---

## 8. What is proved, what is certified, and the named gap

**Proved in full (symbolic / finite-and-complete, no `m`-cap):**
- Lemma 1 closed form; Prop 2 Kummer formula; structural decomposition `Q_4=(a-2)(b-3)H+P_8`;
  the peeling identity.
- `Δ(1) > 0`, `Δ(2) ≥ 0` with the exact `(mod 4)` tie classification, hence `0 ∈ J*` and `|J*|≤2`
  among `j ≤ 7`.
- `Δ(j) > 0` for `3 ≤ j ≤ 7` (five complete residue-class divisibility checks).
- `Δ(j) > 0` for `j ≥ 8` in the **tip-dominated** regime (Case A), via absorption: `Δ(j) ≥
  j-4-2s₂(j-8) ≥ 3`.

**Certified computationally (the one named residual):**
- **(Residual)** `Δ̂(j) > 0` (the heavy-free inequality) for `j ≥ 8` in the **heavy-dominated**
  regime (Case B). Certified: `min Δ̂ = 4,5,6,9,8` for `j ≤ 12`; full closed-form sweep `0`
  violations for `a < 300`. Closing it uniformly needs a `c=4` Number Lemma — a `v₂H(j)` lower
  bound for the sextic heavy quotient `H` (the natural next PROVE target, mirroring how `c=2`'s
  Number Lemma closed its interior).

**Net.** The `c=4` interior tie set is `J* ⊆ {0,2}` (no generator `4/6/8`); proved for all
`1 ≤ j ≤ 7` and `j ≥ 8` (tip-dominated), with the heavy-dominated deep tail certified. This corrects
the earlier overclaim and identifies the precise, structurally understood obstruction to a fully
uniform `c=4` interior theorem.

---

## 9. Diagnosis — why `c=4` is *simpler* than `c=2` in the interior

`c=2`: tip `+P_4(j)=+24\,C(j,4)` (roots `\{0,1,2,3\}`, first nonzero `P_4(4)=24`). The `+24` at the
fourth root `j=4` is what tipped the correction negative there, opening the tie `{0,4}` — a real
generator `4`. `c=4`: tip `+P_8(j)=+8!\,C(j,8)` (roots `\{0,\dots,7\}`). At `j=4` and `j=6` the tip
**vanishes**, so those indices are purely `(a-2)(b-3)H(j)` — no inhomogeneity, no generator. The
first nonzero tip is `P_8(8)=8!`, and the §5 absorption shows even `j=8` stays strictly positive.
So the would-be generators `4,6,8` are all suppressed, and the interior collapses to the single
generator `2`. The general lesson, now visible across `c=1,2,3,4`: `c` sets the tip order
`(2c)!\,C(j,2c)` and hence the *menu* of generators `2,4,\dots,2c`, but which generators actually
fire is a delicate 2-adic event — and for `c=4` only `2` does.

### Files
`projects/code/threerow-c4/`: `c4factor.py` (Lemma 1 derivation), `c4verify.py` (Lemma 1 vs MN +
decomposition), `c4prop2.py` (Prop 2 + wide Δ sweep), `c4peel.py` (peeling identity),
`c4finite.py` (the five §4 residue-class proofs), `c4tail.py` (Case-A bound + wide certification),
`c4caseB.py` (heavy-free `Δ̂`), `c4tie2.py` (`j=2` tie classes), `c4cases.py`/`c4smallj.py`/
`c4comp.py`/`c4Rsign.py`/`c4Hstruct.py`/`c4spec.py` (exploration). MN engine: `threerow-c3/mn.py`.
