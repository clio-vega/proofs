# FINDINGS — Job B: the even-|J*| fpf involution is the TOP-generator toggle, paired by residue mod π²

*Code session 2026-06-14. Script: `2026-06-14-jobB-residue-involution.py`. Engine: `job_jstar_engine.py`.*

## Verdict (plain)

**CONFIRMED on every tie tested (1043 even-|J\*| shapes, m ≤ 16 three-row + all λ⊢2m for m ≤ 11):**

- The XOR-with-top-generator involution on the J\* box is **fixed-point-free** — 0 failures / 1043.
- The leading-π layer `C = Σ_{j∈J*} u_j` cancels **mod π** (the evenness) — 0 failures / 1043.
- On every **|J\*| = 4** tie (124 shapes), the **top-generator pairing** `{j, j+4}` is the **unique**
  pairing whose partners cancel **mod π²** (`v_π(sum) ≥ 2`); Sawin's adjacent `j↔j+2` and the
  outer/cross pairing each leave `v_π(sum) = 1` — **0 / 124** for adjacency. Sawin adjacency is
  **FALSIFIED** as the even-|J\*| involution; the top-generator toggle is confirmed.

No counterexample up to m = 16. (One apparent anomaly at `(24,5,3)` was a `sympy.nsimplify`
artifact corrupting the Gaussian integer `u_0 = −6513255 i`; the raw `expand(term/π^V)` is already a
clean Gaussian unit. Fixed — see "honesty note" below.)

## Setup (the verified identity that makes the test faithful)

With `π = 1+i` and `M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩`:

```
G_λ(i) = ⟨s_λ,(p₂+π e₂)^m⟩ = Σ_j C(m,j) M_j (i−1)^j ,   v_π(term_j) = val(j) = j + 2 v₂(C(m,j)M_j).
```
`J* = argmin_j val(j)`, `V = min val`. The **leading unit** is `u_j = C(m,j)M_j(i−1)^j / π^V`. Because
`i−1 = i·π` exactly, `u_j = C(m,j)M_j · i^j / π^{V−j}` is always **±odd or ±odd·i** — a genuine
Gaussian unit (`v_π = 0`).

## The residue-mod-π² law (the mechanism, sharper than the note had it)

`π² = (1+i)² = 2i`, and the ideal `(2i) = (2)`, so `ℤ[i]/(π²) = ℤ[i]/(2)` has 4 classes and
`residue(u) = (Re u mod 2) + (Im u mod 2)·i`. A unit has **odd norm**, so exactly one of Re, Im is
odd ⟹ `residue(u) ∈ {1, i}`.

**Law (0 failures / 1043):** within any J\*, `residue(u_j) = residue(u_{j'})  ⟺  j ≡ j' (mod 4)`.

So residue mod π² is governed by `j mod 4` (up to a global per-shape swap of `{1,i}` — e.g. `(9,6,3)`
has `j≡3 mod4 → 1`, `(17,10,3)` has `j≡3 mod4 → i`; the *relative* residues are fixed by `j mod 4`).
Consequences:

- **shift by 4** (≡ 0 mod 4) → residue-**equal** → `u + u = 2r`, `v_π(2r) = 2` → **cancels mod π²**.
- **shift by 2** (Sawin adjacency) → residue-**opposite** `(1, i)` → `u + u ~ 1+i = π`, `v_π = 1` →
  **does not** cancel mod π².

This is *why* the top generator governs: the J\* box is `j₀ + span_{GF(2)}{2^a : a∈S}` (generators
2, 4, …), and XOR-with-the-top-generator shifts by the largest `2^a`. For |J\*| = 4 (`S = {1,2}`,
box `{3,5,7,9}=3+{0,2,4,6}`) the top shift is **4** → residue-equal. For Sawin's smallest-generator
flip the shift is **2** → residue-opposite → no mod-π² cancellation.

## Representative |J*| = 4 residue table (the discriminating case)

```
λ=(9,6,3)  m=9  V=17  J*={3,5,7,9}
   j=3: u=21021       res=1     j=7: u=-549     res=1
   j=5: u=-8001 i     res=i     j=9: u=3 i      res=i
   TOP-gen {3,7},{5,9}:  v_π(sum)=6, 2    -> CANCELS (residue-equal)   <== the involution
   Sawin   {3,5},{7,9}:  v_π(sum)=1, 1    -> no cancel (residue 1 vs i)
   cross   {3,9},{5,7}:  v_π(sum)=1, 1    -> no cancel
```
Same pattern on all 8 sampled |J\*|=4 shapes incl. `(13,6,3),(13,10,3),(17,6,3),(17,10,3),(21,6,3)`
(three-row c=3 family) and non-three-row `(6,3,1,1,1),(5,2,2,1,1,1)`.

## The |J*| = 2 nuance (honest caveat — needed for the involution program)

The basic evenness needs only **mod π**: every unit ≡ 1 ≡ i (mod π) because `1−i = −i·π`, so **any**
pairing sums to 0 mod π once |J\*| is even. The **mod π²** refinement only has content when |J\*| ≥ 4.
For |J\*| = 2 the single pair reaches mod π² **iff** its box generator is 4 (residue-equal); the common
gap-2 box (generator 2) cancels only mod π (`v_π(C) = 1`). This does **not** weaken the involution
claim — for |J\*| = 2 there is only one pairing and "which pairing" is vacuous — but it sharpens the
spec for the Pfaffian home: the residue-class pairing it must realize is **shift-by-a-multiple-of-4**
(top generator), and it is only at |J\*| ≥ 4 that this is a non-trivial constraint distinguishing it
from nearest-neighbour matchings.

## Implication for the surviving Pfaffian home (Rains–Warnaar / Fischer–Gangl at t=i)

The involution's currency is **leading-unit residue mod π² ≈ j mod 4**, and its action is XOR with the
**top** generator (shift by a multiple of 4), *not* nearest-neighbour adjacency. Any Pfaffian
realization specialised at `t = i` must produce surviving sectors that pair the **residue-equal**
indices `{3,7},{5,9}` (shift 4), not the adjacent `{3,5},{7,9}` (shift 2). This is the sharpened
target (optional Pfaffian-sector check was not run — no readily-computable bounded-Littlewood data
in-container without Sage; reported as skip).

## Honesty note

`sympy.nsimplify` must **not** be applied to `expand(term/π^V)`: on `(24,5,3)` it turned the exact
Gaussian unit `−6513255 i` into an irrational expression, producing a spurious `residue 0` and a
false "mod-π cancellation failure." The exact quotient is already a Gaussian integer; the fix casts
`re/im` to `int` and asserts integrality. After the fix: **0 failures** across all four checks.
