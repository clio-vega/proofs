# FINDINGS — Job B: three-row M_j + sharp J* box (2026-06-13 code)

**Script:** `2026-06-13-jobB-threerow.py` (reuses `job_jstar_engine.py`).
Scope: all three-row `λ=(a,b,c)`, `a≥b≥c≥1`, `λ⊢2m`, `m ≤ 12` (2226 `(λ,j)` pairs, 214 shapes).
`M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩`, `val(j) = j + 2v₂(C(m,j)M_j)`, `J* = argmin`.

## Result 1 — the exact closed form: **signed S₃ JT determinant of trinomials** ✓

`M_j = Σ_{w∈S₃} sgn(w) · T_j(λ+δ − w·δ)`,  `δ=(2,1,0)`,
where `T_j(α) = [s^{α₁}t^{α₂}u^{α₃}] e₁^{2(m−j)} e₂^j` (the trinomial expansion of `e₁²,e₂` in
3 variables). This is the GL₃ Weyl-character extraction of the `s_λ`-coefficient of
`(e₁²+x e₂)^m`, since `Φ_λ(x)=⟨s_λ,(e₁²+xe₂)^m⟩ = Σ_j C(m,j) M_j x^j`.

**Cross-check: 0 failures / 2226 pairs** against the vertical-2-strip chain engine.
The trinomial determinant **is** `M_j` — this is the load-bearing closed form for PROVE.

## Result 2 — the naive dimension conjecture is **FALSE** for three rows ✗

`M_j ?= f^{(a−j,b−j,c−j)}` (3-row dim, all rows shortened by `j`):
- **HOLD 634 / 2226, BREAK 1592 / 2226.**
- Holds at `j=0` always (`M_0 = f^λ`, the only universally safe case).
- **Breaks even when `c−j≥0`** (the shift is a genuine partition): hold 214 / break 560.
- Smallest break: `λ=(2,1,1), j=1`: `M_1 = 2` but `f^{(1)} = 1`. Then `(2,2,2)`: `M_1=3 vs 1`;
  `(3,2,1)`: `M_1=8 vs 2`, `M_2=4 vs 0`.

So the clean two-row identity `M_j = f^{(a−j,b−j)}` (memory `2026-06-12-tworow-Jstar-even`) does
**NOT** lift to three rows. PROVE must use the **signed determinant**, not a single dimension.
(The vertical-2-strip removal `λ→μ` lands in many shapes `μ⊢2(m−j)`, not just the row-shift; the
2-row case collapses because there is only one chain, three rows do not.)

## Result 3 — the sharp three-row J* box: **{0,2,4,6}**, generators {2,4}

`|J*|` census (214 three-row shapes, `m ≤ 12`):

| \|J*\| | 1 | 2 | 4 |
|---|---|---|---|
| count | 148 | 64 | 2 |

- **All `J*` single-parity** (0 mixed) → half-poly split always valid.
- **All offset sets `J*−minJ*` are XOR-closed 2-adic boxes** (0 failures).
- Distinct offset patterns: `{0}`×148, `{0,2}`×42, **`{0,4}`×22**, **`{0,2,4,6}`×2**.
- So the sharp three-row box is `J*−min ⊆ {0,2,4,6} = 2·{0,1,2,3} = ⟨2,4⟩_{XOR}`, generators
  among `{2,4}`; `|J*| = 2^{|S|}`, `S⊆{2,4}`.
- **The strong box `{0,2}` is FALSE past `c=1`.** The new generator `4` is real: `{0,4}` already
  appears at `|J*|=2` (22 shapes), and the full box `{0,2,4,6}` at `|J*|=4`.

**First `|J*|=4` among three-row shapes:** `m=9, λ=(9,6,3), J*={3,5,7,9}` (offsets `{0,2,4,6}`).
Second: `m=11, λ=(13,6,3)`. Note `(9,6,3)=3·(3,2,1)` — a scaled staircase; the `J*={3,5,7,9}`
is **odd-parity** (`a=9` odd), matching the parity rule of `2026-06-13-threerow-c1-Jstar-even`
(`a` odd → odd-parity `J*`) now seen to persist at `c=3`.

## PROVE targets pinned by this run

1. The right object is the **S₃ trinomial determinant**, not a shifted dimension. The 2-row
   dimension shortcut does not generalise.
2. The correct sharp box to prove for three rows is **`{0,2,4,6}` (generators 2 and 4)**, not
   `{0,2}`. Proving even-`|J*|` for three rows means proving the offset set is an XOR box in
   `2·{0,1,2,3}`.
3. `|J*|=4` debut at the scaled staircase `(9,6,3)` is the minimal witness — the place to test any
   candidate generator-`4` toggle / fpf-involution argument.
