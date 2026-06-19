# The c=4 interior Number Lemma: `16 ∣ H` closes the §6 residual ⟹ three-row `(a,b,4)` is COMPLETE

**Date:** 2026-06-19 (prove session)
**Status:** The last named residual of the three-row `(a,b,4)` family is **closed**. Together with the
already-complete `c=4` boundary (`2026-06-16-c4-boundary-complete.md`) and the interior indices
`1 ≤ j ≤ 7` (`2026-06-18-c4-interior-Jstar-even.md`, §§3–4), this makes `(a,b,4)` a **genuinely
complete** family — interior *and* boundary — joining `c=1,2,3`.

This note proves the deep-index inequality `Δ(j) > 0` for all interior `j ≥ 8`, which was previously
only *certified* (`m ≤ 120`, `a < 300`). The engine is a single clean lemma:

> **Lemma N4 (c=4 Number Lemma).** `16 ∣ H(a,b,j)` for all integers `a ≡ b (mod 2)` and all `j`.
> Equivalently `v₂H(j) ≥ 4`. (Sharp: `32 ∤ H` for some residues.)

`H` is the sextic heavy quotient of the structural decomposition
`Q_4(j) = (a−2)(b−3)H(j) + P_8(j)`. The §6 residual asked for exactly such a `v₂H` floor; the
surprise is that the floor is the **constant `4`**, provable by a finite residue check, and that this
constant — fed through the §5 absorption — suffices for every `j ≥ 8` save two indices `j ∈ {8,10}`,
which fall to direct finite residue checks of the same kind that closed `3 ≤ j ≤ 7`.

All notation is that of `2026-06-18-c4-interior-Jstar-even.md`. Recall (proved there):

- `Δ(j) = j − 2s₂(j) + 2T(j)`, `T(j) = v₂C(a+5,j) + v₂C(b+4,j) + v₂Q_4(j) − v₂Q_4(0)` (Prop 2);
- `Q_4(j) = (a−2)(b−3)H(j) + P_8(j)`, `P_8(j) = (j)_8 = 8!·C(j,8)`, `H` sextic in `j`,
  `H(0) = (a+3)(a+4)(a+5)(b+2)(b+3)(b+4)`;
- **Peeling identity** (for `j ≥ 3`):
  `T(j) = v₂((a+2)_{j−3}) + v₂((b+1)_{j−3}) − 2·vfact(j) + v₂Q_4(j) − v₂(a−2) − v₂(b−3)`,
  where `(x)_r = x(x−1)···(x−r+1)` and `vfact(r) = v₂(r!) = r − s₂(r)`.

Throughout `a ≥ b ≥ 4`, `a ≡ b (mod 2)`, and `j` is interior: `8 ≤ j ≤ b`. Since `a ≡ b (mod 2)`,
exactly one of `a−2, b−3` is even; set `φ := v₂((a−2)(b−3)) = v₂(a−2) + v₂(b−3)`. As `a ≥ 4`, the
even one is `≥ 2`, so **`φ ≥ 1`**.

---

## 1. Lemma N4: `16 ∣ H`

**Statement.** For all integers `a ≡ b (mod 2)` and all integers `j`, `H(a,b,j) ≡ 0 (mod 16)`.

**Proof.** `H` is a polynomial with integer coefficients, so for any modulus `N` the residue
`H(a,b,j) mod N` depends only on `(a mod N, b mod N, j mod N)` (each monomial satisfies
`(x+N)^k ≡ x^k (mod N)`). Hence `H mod 16` is a function of `(a,b,j) mod 16`, and the congruence
`16 ∣ H` for all `a ≡ b (mod 2)` is equivalent to the **finite** assertion that
`H(a,b,j) ≡ 0 (mod 16)` for every residue triple `a,b,j ∈ {0,…,15}` with `a ≡ b (mod 2)` (the
parity constraint is preserved mod 16). The exhaustive check (`c4N4.py`, all `16·8·16 = 2048`
admissible triples) returns **0 exceptions**. Therefore `v₂H(j) ≥ 4`. ∎

*Sharpness.* `32 ∤ H` for some admissible residue triple, so the floor `4` cannot be raised:
e.g. `v₂H(10,8,8) = 4` (and `v₂H(7,7,0) = 4`); the minimum `v₂H = 4` is attained for every `j` in the
tested range. The bound is exactly the constant `4`, uniform in `a,b,j`.

A conceptual remark. `H(j) − H(0) = −2j·R(a,b,j)` with `R ∈ ℤ[a,b,j]` (so `H(j) ≡ H(0) (mod 2)`,
and as `H(0)` is a product of three consecutive integers times three consecutive integers it is
always even — the crude `v₂H ≥ 1`). Lemma N4 sharpens this all the way to `16`; unlike the `c=2`
quartic, no `a,b`-dependent factor of `H` is needed — `H` carries an absolute 2-adic floor.

---

## 2. The case split on `v₂Q_4(j)`

Fix interior `j ≥ 8`. Compare the two summands of `Q_4(j) = (a−2)(b−3)H(j) + P_8(j)`. Write
`P := v₂(P_8(j))` and `Hv := v₂((a−2)(b−3)H(j)) = φ + v₂H(j)`. By Kummer
`v₂C(j,8) = s₂(8) + s₂(j−8) − s₂(j) = 1 + s₂(j−8) − s₂(j)`, so

> `P = v₂(8!) + v₂C(j,8) = 7 + (1 + s₂(j−8) − s₂(j)) = 8 + s₂(j−8) − s₂(j)`,  and  `P ≥ 7`.

### Case A: `Hv ≥ P`  (heavy dominates, `v₂Q_4(j) ≥ P`)

This is the §5 absorption argument; we reproduce it for completeness. For `j ≥ 8` the run
`(a+2)_{j−3}` (length `j−3 ≥ 5`, all factors positive since `a ≥ b ≥ j ≥ 8` ⟹ bottom factor
`a−j+6 ≥ 6`) splits as

`  (a+2)_{j−3} = [(a+2)(a+1)a(a−1)]·(a−2)·(a−3)_{j−8}`,

a 4-block (`v₂ ≥ vfact(4) = 3`), the single factor `a−2` (`v₂ = v₂(a−2)`), and a `(j−8)`-block
(`v₂ ≥ vfact(j−8)`). Hence `v₂((a+2)_{j−3}) ≥ 3 + v₂(a−2) + vfact(j−8)`, and symmetrically
`v₂((b+1)_{j−3}) ≥ 3 + v₂(b−3) + vfact(j−8)` (bottom factor `b−j+5 ≥ 5`). Substituting in the
peeling identity (the `v₂(a−2), v₂(b−3)` cancel):

`  T(j) ≥ 6 + 2·vfact(j−8) − 2·vfact(j) + v₂Q_4(j) = 6 − 2P + v₂Q_4(j) ≥ 6 − P`

(using `vfact(j) − vfact(j−8) = v₂P_8(j) = P` and `v₂Q_4(j) ≥ P` in Case A). Therefore

`  Δ(j) = j − 2s₂(j) + 2T(j) ≥ j − 2s₂(j) + 12 − 2P = j − 2s₂(j) + 12 − 2(8 + s₂(j−8) − s₂(j))`
`        = j − 4 − 2s₂(j−8).`

Writing `j = 8 + t` (`t ≥ 0`): `Δ(j) ≥ 4 + (t − 2s₂(t)) ≥ 3`, since `2s₂(t) ≤ t + 1` for all
`t ≥ 0` (Lemma A; equality at `t = 1,3`). **So `Δ(j) ≥ 3 > 0`.** ∎ (Case A)

### Case B: `Hv < P`  (tip strictly larger valuation)

Then `v₂Q_4(j) = Hv = φ + v₂H(j)` exactly, and the peeling identity collapses (the `−v₂(a−2)−v₂(b−3)
= −φ` cancels the `+φ`):

> `  Δ(j) = Δ̂(j) := j − 2s₂(j) + 2[ v₂((a+2)_{j−3}) + v₂((b+1)_{j−3}) − 2·vfact(j) + v₂H(j) ].`

This is precisely the heavy-free quantity of §6. We bound it two ways.

**Sub-case B1: `j ∉ {8,10}`.** Use the elementary `v₂((x)_r) ≥ vfact(r)` (a product of `r`
consecutive integers is divisible by `r!`) on both falling factorials, and Lemma N4 `v₂H(j) ≥ 4`:

`  Δ̂(j) ≥ j − 2s₂(j) + 2[ 2·vfact(j−3) − 2·vfact(j) + 4 ] = j − 2s₂(j) − 4(vfact(j) − vfact(j−3)) + 8.`

Since `vfact(j) − vfact(j−3) = v₂(j(j−1)(j−2))`, this is

> `  Δ̂(j) ≥ g(j) := j − 2s₂(j) − 4·v₂(j(j−1)(j−2)) + 8.`

**Claim: `g(j) > 0` for every `j ≥ 8` with `j ∉ {8,10}`.**

*Bound on the valuation.* Among `j, j−1, j−2` at most one contributes a high power of `2`:
- `j` odd: only `j−1` even, `v₂(j(j−1)(j−2)) = v₂(j−1) ≤ ⌊log₂ j⌋`.
- `j` even: `j−1` odd; `j, j−2` are consecutive even numbers, so exactly one is `≡ 2 (mod 4)`
  (valuation `1`), whence `v₂(j) + v₂(j−2) = 1 + max(v₂(j), v₂(j−2)) ≤ 1 + ⌊log₂ j⌋`.

In both cases `v₂(j(j−1)(j−2)) ≤ 1 + ⌊log₂ j⌋`. With `s₂(j) ≤ ⌊log₂ j⌋ + 1`,

`  g(j) ≥ j − 2(⌊log₂ j⌋ + 1) − 4(1 + ⌊log₂ j⌋) + 8 = j − 6⌊log₂ j⌋ + 2 ≥ j − 6log₂ j + 2.`

The real function `h(x) = x − 6log₂ x + 2` has `h'(x) = 1 − 6/(x ln 2) > 0` for `x > 6/ln2 ≈ 8.66`,
and `h(32) = 32 − 30 + 2 = 4 > 0`; hence `g(j) > 0` for all `j ≥ 32`. For the finite range
`11 ≤ j ≤ 31` the values are computed directly (`c4caseB1.py`):

`  g(11..20) = 9, 4, 7, 4, 11, 2, 5, 2, 17, 12,  …`  (all `> 0`),

and the only `j ≥ 8` with `g(j) ≤ 0` are `j = 8` (`g = −2`) and `j = 10` (`g = −2`). This proves the
claim, so **`Δ̂(j) ≥ g(j) > 0`** for `j ∉ {8,10}`. ∎ (Case B1)

**Sub-case B2: `j ∈ {8,10}`.** Write `Π_j := (a+2)_{j−3}·(b+1)_{j−3}·H(j)`, so
`Δ̂(j) = j − 2s₂(j) + 2[v₂(Π_j) − 2·vfact(j)]`. Both `j` are even, so `Δ̂(j)` is even, and `Δ̂(j) > 0`
iff `Δ̂(j) ≥ 2`:

- `j = 8`: `vfact(8) = 7`, `Δ̂(8) = 6 + 2v₂(Π_8) − 28 = 2v₂(Π_8) − 22`. Positive iff `v₂(Π_8) ≥ 12`,
  i.e. **`2¹² ∣ Π_8`**.
- `j = 10`: `vfact(10) = 8`, `Δ̂(10) = 6 + 2v₂(Π_{10}) − 32 = 2v₂(Π_{10}) − 26`. Positive iff
  `v₂(Π_{10}) ≥ 14`, i.e. **`2¹⁴ ∣ Π_{10}`**.

Each is a polynomial divisibility in `(a,b)` over the lattice `a ≡ b (mod 2)`, hence a finite residue
check modulo `2¹²` resp. `2¹⁴`. The exhaustive checks (`c4N4.py`):

- `2¹² ∣ (a+2)_5 (b+1)_5 H(8)` for all `a,b mod 2¹²` with `a ≡ b (mod 2)`: **0 failures**;
- `2¹⁴ ∣ (a+2)_7 (b+1)_7 H(10)` for all `a,b mod 2¹⁴` with `a ≡ b (mod 2)`: **0 failures**.

(These are complete residue-class proofs of the same type as the §4 checks for `3 ≤ j ≤ 7`; the
larger moduli are exactly `2^{k_j}` with `k_8 = 12, k_{10} = 14`.) Therefore `Δ̂(j) ≥ 2 > 0`. ∎ (B2)

---

## 3. Conclusion

Combining: for every interior `j ≥ 8` either Case A (`Δ(j) ≥ 3`) or Case B holds, and in Case B
`Δ(j) = Δ̂(j) > 0` by B1 (`j ∉ {8,10}`) or B2 (`j ∈ {8,10}`). Hence

> **Theorem.** For every three-row shape `λ = (a,b,4) ⊢ 2m` and every interior index `8 ≤ j ≤ b`,
> `Δ(j) > 0`.

Together with `2026-06-18-c4-interior-Jstar-even.md` §§3–4 (`Δ(1) > 0`, `Δ(2) ≥ 0` with the exact
`(mod 4)` tie classification, `Δ(j) > 0` for `3 ≤ j ≤ 7`), this establishes the full interior:

> `0 ∈ J*`, `J* ⊆ {0,2}`, with `J* = {0,2}` iff `(a,b) ≡ (2,2)` or `(1,3) (mod 4)`. In every tie
> `|J*| ≤ 2` is even — i.e. `G_λ(i) = 0 ⟹ |J*|` even for the `c=4` interior, **with no residual.**

Adjoining the boundary lemma (`2026-06-16-c4-boundary-complete.md`, already complete), the three-row
family **`(a,b,4)` is now fully proven, interior and boundary**, on the same footing as `c=1,2,3`.

---

## 4. What is proved (no residual remains)

**Proved in full (symbolic / complete finite residue checks, no `m`-cap):**
- Lemma N4: `16 ∣ H` (residues mod 16) — sharp.
- Case A (`v₂Q_4 ≥ v₂P_8`): `Δ(j) ≥ j − 4 − 2s₂(j−8) ≥ 3` (§5 absorption, falling-factorial split).
- Case B1 (`j ∉ {8,10}`): `Δ̂(j) ≥ g(j) > 0`, with `g(j) > 0` proved for `j ≥ 32` analytically and
  `11 ≤ j ≤ 31` by direct evaluation; `g(j) ≤ 0` only at `j ∈ {8,10}`.
- Case B2 (`j ∈ {8,10}`): `2¹² ∣ Π_8`, `2¹⁴ ∣ Π_{10}` (complete residue checks) ⟹ `Δ̂(j) ≥ 2`.

**Net.** The §6 residual ("heavy-free `Δ̂(j) > 0` for `j ≥ 8`, certified not proved") is closed. The
key was Lemma N4 — the `v₂H` floor the §6 note called for — turning out to be the *constant* `4`,
which (with the absorption identity) clears all but two indices, themselves finite checks.

### Cross-validation (against the closed form / MN ground truth)
`c4crossval.py`: for all `λ=(a,b,4)`, `a ≤ 160`, every interior `j ≥ 8` (301 378 triples): each
triple is classified Case A (132 638) / Case B (168 740), `Δ(j) > 0`, and each lies within the bound
its case proves (Case A `≥ 3`; B1 `≥ g(j)`; B2 `≥ 2`). **0 problems.** Moreover `min(Δ − g) = 0`
over Case B1 — the crude bound is attained, so `g` is the exact floor (not loosenable to drop the
`{8,10}` checks). Lemma N4 sharp (`v₂H(10,8,8) = 4`, floor `4` attained for every `j`). The `Δ̂`
minima per `j` reproduce
`4,5,6,9,8,…` exactly.

### Files
`projects/code/threerow-c4/`: `c4N4.py` (Lemma N4 mod-16 proof + the `j=8` `2¹²` and `j=10` `2¹⁴`
residue checks), `c4caseB1.py` (the `g(j)` inequality + threshold), `c4crossval.py` (end-to-end case
split vs closed form). Reuses `H` from `c4finite.py`, MN engine `threerow-c3/mn.py`.
