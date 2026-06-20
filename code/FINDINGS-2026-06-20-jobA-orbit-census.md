# FINDINGS — Job A (2026-06-20 CODE): orbit-size census of the 2-adic toggle  [for Rick]

**Script:** `jobA_orbit_census.py` · **Ground truth:** Murnaghan–Nakayama `M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩`
(abacus χ engine, `threerow-c3/mn.py`). Every `|J*|` is computed **twice** — once from MN, once from
the proved closed form — and the two are asserted equal.

## One-line answer for Rick

The toggle's orbit is **all-or-nothing on the residue class**: inside every `(a,b) mod 4` class where
the gate *permits* a tie, **100 %** of shapes tie (`|J*|=2`); in every other class **0 %** do (`|J*|=1`).
This is the abundance analogue of your 17/17, 21/21 — here it is **`364/364` (c=4)** and **`325/325` (c=5)**,
exact, over the interior range `m ≤ 60`, `a > b > c`.

- **c=4:** 1431 interior shapes, **364 toggles / 1067 rigid**, MN-vs-closed mismatches **0**, box violations **0**.
- **c=5:** 1378 interior shapes, **325 toggles / 1053 rigid**, MN-vs-closed mismatches **0**, box violations **0**.
- Single-generator confirmed: every `J*` is either `{j₀}` or `{j₀, j₀+2}` — never larger, step always 2.

## c=4 census table  `(a%4, b%4) → orbit class`  (`m ≤ 60`)

```
 (a%4,b%4) | #|J*|=1 | #|J*|=2 | j₀ | tie-gate
   (0,0)   |   169   |     0   |  0 |
   (0,2)   |   182   |     0   |  0 |
   (1,1)   |   182   |     0   |  0 |
   (1,3)   |     0   |   182   |  0 |  TIE
   (2,0)   |   169   |     0   |  0 |
   (2,2)   |     0   |   182   |  0 |  TIE
   (3,1)   |   196   |     0   |  0 |
   (3,3)   |   169   |     0   |  0 |
```
- **j₀ = 0 for every class** (c even ⟹ no forced descent). Tie classes exactly `{(2,2),(1,3)}`.
- Minimal tie witnesses: `(a,b,m)=(10,6,10)` `J*={0,2}` for (2,2); `(9,7,10)` `J*={0,2}` for (1,3).

## c=5 census table  `(a%4, b%4) → orbit class`  (`m ≤ 60`)

```
 (a%4,b%4) | #|J*|=1 | #|J*|=2 | j₀ | tie-gate
   (0,1)   |     0   |   156   |  0 |  TIE
   (0,3)   |   182   |     0   |  0 |
   (1,0)   |   169   |     0   |  3 |
   (1,2)   |   182   |     0   |  3 |
   (2,1)   |   169   |     0   |  0 |
   (2,3)   |   169   |     0   |  0 |
   (3,0)   |     0   |   169   |  3 |  TIE
   (3,2)   |   182   |     0   |  3 |
```
- **The offset is set by parity(a):** `a even ⟹ j₀ = 0`, `a odd ⟹ j₀ = 3` (the forced descent −1,−2,−3
  for c odd). Confirmed in the `j₀` column: every odd-`a` class shows `j₀=3`, every even-`a` class `j₀=0`.
- Tie classes exactly `{(0,1) [a even, box {0,2}], (3,0) [a odd, box {3,5}]}`.
- Minimal tie witnesses: `(a,b,m)=(12,9,13)` `J*={0,2}` for (0,1); `(11,8,12)` `J*={3,5}` for (3,0).

## The object to cite

| c | orbit size `|J*|` | shapes | mod-4 class(es) | tie fraction in gate |
|---|---|---|---|---|
| 4 | 1 (rigid) | 1067 | (0,0),(0,2),(1,1),(2,0),(3,1),(3,3) | — |
| 4 | 2 (toggle) | 364 | (1,3),(2,2) | **364/364 = 1** |
| 5 | 1 (rigid) | 1053 | (0,3),(1,0),(1,2),(2,1),(2,3),(3,2) | — |
| 5 | 2 (toggle) | 325 | (0,1),(3,0) | **325/325 = 1** |

**Verdict.** The toggle is a *rigid* residue-class indicator, not a fractional one — exactly the c=3
behaviour, and the single-generator structure (`|J*| ≤ 2`) of c=4/c=5 makes the per-shape orbit trivially
1 or 2 with **no shape in a tie class failing to tie**. Cross-check: **0** MN-vs-closed-form mismatches over
all 2809 interior shapes.
