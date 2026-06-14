# FINDINGS — Job B: c=3 closed form, J* box census, Sawin involution test

**Date:** 2026-06-14 (system clock) · code session · feeds the parallel c=3 PROVE
**Scripts:** `2026-06-14-jobB-census-sawin.py` (main), `…-jobB-fit.py`, `…-jobB-c3-explore.py`
**Engine:** `job_jstar_engine.py` — `M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩`, `val(j)=j+2v₂(C(m,j)M_j)`, `J*=argmin`.

---

## Summary of verdicts

| Item | Result |
|---|---|
| **1. c=3 closed form** | c=2 form **verified** (5883 checks, 0 fail). The single-binomial × sextic collapse **DOES hold for c=3** — the PROVE session's `c3factor.py` Lemma 1, which I **independently verified against my engine (1627 checks, 0 fail)**: same binomial `C(2(m−j),b−j)`, prefactor `(a−b+1)`, sextic `Q₃`, denominator `6(a+4−j)(b+1−j)(b+2−j)(b+3−j)`. (My first attempt wrongly reused c=2's denominator `(a+3−j)(b+2−j)(b+1−j)` → degree grew with b; the c=3 denominator has **a+4 and three b-factors b+1,b+2,b+3**.) |
| **2. Generator-S congruence** | `|J*|∈{1,2,4}`, `S⊆{gen 2, gen 4}`, and **S is determined by (a,b) mod 4** (0 ambiguous classes). **The modulus does NOT grow from c=2.** |
| **3. Sawin involution** | **NO.** Consecutive (= gen-2, `j↔j+2`) pairing does not cancel — adjacent leading units have opposite residues `1` vs `i`. The actual fpf involution is the **top-generator toggle `j↔j+4`** (equal-residue pairs, `r+r=2r≡0 mod π²`). |

---

## Item 1 — c=2 closed form verified; c=3 closed form independently confirmed

**c=2 closed form (binomial has a j-dependent top, like the two-row case):**
```
M_j = (a−b+1) · C(2(m−j), b−j) · Q(j) / [ 2 (a+3−j)(b+2−j)(b+1−j) ],   valid for 0 ≤ j ≤ b
Q(j) = a(b−1)[(a+3)(b+2) − 2j²] + j(j−1)(j−2)(j−3)        (quartic; inhomog. part = falling fact. j^{(4)})
```
**Verified: 5883 (a,b,2)-and-j checks, 0 failures.** (My earlier failed search wrongly fixed `N`; the
top is `2(m−j)=a+b+2−2j`, not constant.)

**c=3 closed form (the PROVE session's Lemma 1 in `c3factor.py`), independently verified here:**
```
M_j = (a−b+1) · C(2(m−j), b−j) · Q₃(j) / [ 6 (a+4−j)(b+1−j)(b+2−j)(b+3−j) ],   N=2(m−j)=a+b+3−2j
Q₃ = −j⁶ + 15j⁵ − 79j⁴ + 201j³ − 244j² + 36j + 144      (pure-j part; full Q₃ has a,b-coefficients)
```
**Independently verified against my engine: 1627 checks (all (a,b,3), all j, m≤20), 0 failures** — a
genuine cross-check via a different code path (my `M_j` chain vs their `D`-sum derivation). The single-
binomial × sextic collapse **does hold for c=3**; the leading `j⁶` coefficient is `−1`, and structurally
`Q₃ = (a−1)(b−2)·H − 720 C(j,6)` (their `c3decomp.py`), with `(a−1)(b−2)=(a−c+2)(b−c+1)|_{c=3}` and the
inhomogeneous falling factorial `720 C(j,6) = j^{(6)}` (6 roots {0,…,5}) — the exact c-pattern lift of
c=1's `j^{(2)}` and c=2's `j^{(4)}`.

**Lesson / why my first attempt missed it:** I reused c=2's denominator `(a+3−j)(b+2−j)(b+1−j)`; the
correct c=3 denominator is `6·(a+4−j)·(b+1−j)(b+2−j)(b+3−j)` — `a+4` (not `a+3`) and **three** consecutive
b-factors (not two). With the wrong denominator the residual is non-polynomial (degree → b); with the right
one it is exactly the sextic. (Note: like c=2, the form is the leading-binomial representation for `j≤b`;
`v_π`-relevant `J*` can sit at `j>b` — e.g. `(9,6,3)`: `J*={3,5,7,9}`, `b=6` — where the determinant /
the analytically continued binomial supplies the value.)

**v₂(Q₃(0)) (the c=3 analogue of c=2's `v₂Q(0)` for Prop-2):** varies — `(9,6,3)→10`, `(8,5,3)→5`,
`(13,10,3)→11`. So the Prop-2 base term `−v₂Q₃(0)` is shape-dependent, as at c=2.

## Item 2 — the c=3 box census (m ≤ 18, all (a,b,3) with a+b odd)

`|J*|` distribution: **{1: 68, 2: 28, 4: 9}**. Generator sets (offsets are `2·{half}`):

| S (generators) | meaning | count | **congruence rule (a,b) mod 4** |
|---|---|---|---|
| `{}` | `|J*|=1` | 68 | (0,3),(1,0),(2,1),(3,0),(3,2) |
| `{gen 2}` → offsets {0,2} | `|J*|=2` | 16 | **(a,b)≡(2,3)** |
| `{gen 4}` → offsets {0,4} | `|J*|=2` | 12 | **(a,b)≡(0,1)** |
| `{gen 2, gen 4}` → {0,2,4,6} | `|J*|=4` | 9 | **(a,b)≡(1,2)** |

**S is a clean function of (a,b) mod 4 — 0 ambiguous classes** (mod 8, mod 16 are consistent refinements;
mod 4 already suffices). **The modulus does not grow from c=2** (also mod 4). For the whole `|J*|=4`
family, **J\* = {3,5,7,9} = 3 + {0,2,4,6} is constant** (independent of a,b). Minimal witness
`(9,6,3)=3·(3,2,1)`. First gen-4-only shape: `(8,5,3)` → `J*={0,4}`.

## Item 3 — Sawin adjacent-pair involution: NO; the involution is j↔j+4

**Derived & verified identity** (makes the leading-π test faithful):
```
G_λ(i) = ⟨s_λ,(p₂+π e₂)^m⟩ = Σ_j C(m,j) M_j (i−1)^j,   π=1+i,
   since p₂ = h₁² − 2e₂ ⟹ p₂+πe₂ = h₁² + (i−1)e₂, and ⟨s_λ,(h₁²+x e₂)^m⟩ = Σ_j C(m,j)M_j x^j.
   v_π(term_j) = j + 2v₂(C(m,j)M_j) = val(j).   J* = argmin = minimal-π-valuation terms.
```
Verified `v_π(G_λ(i)) = min val(J*)` when `|J*|` odd, and `> min val` (leading layer cancels) when `|J*|`
even — e.g. `(9,6,3)`: `G=−43008i`, `v_π(G)=22 > 17=min val`, `|J*|=4`.

**The leading units `u_j = C(m,j)M_j(i−1)^j / π^V` (j∈J\*) and their residues mod π²(=mod 2):**
for every `|J*|=4` shape, residues **alternate** along sorted J\* — `j≡3 (mod 4) → one of {1,i}`,
`j≡1 (mod 4) → the other`. Concretely at `(9,6,3)`: `u=(21021, −8001i, −549, 3i)` for `j=(3,5,7,9)`.

- **Sawin's consecutive (adjacent) pairing `(j,j+2)`** pairs **opposite** residues (`1` with `i`):
  `1+i` has `v_π=1` — does **not** cancel the layer. **Sawin's rule fails here.**
- **The top-generator toggle `j↔j+4`** pairs **equal** residues `{3,7}` and `{5,9}`:
  `r+r = 2r`, `v_π=2` — cancels mod π². **This is the fixed-point-free involution forcing `|J*|` even.**

So the even-`|J*|` pairing for c=3 is the **largest** generator toggle (`+4`), **not** the smallest/adjacent
one the Sawin and earlier "smallest-generator" guesses predicted. Useful redirect for the involution program:
model the **top** generator, and the pairing is by leading-unit residue mod π² (≈ `j mod 4`), not by adjacency.

## Files
- `2026-06-14-jobB-census-sawin.py` — identity verification, box census, congruence rule, Sawin test.
- `2026-06-14-jobB-fit.py`, `…-c3-explore.py`, `…-denomsearch.py` — closed-form fitting (exploratory).
