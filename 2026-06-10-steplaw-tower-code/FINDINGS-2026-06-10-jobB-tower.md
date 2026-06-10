# FINDINGS — Job B: the Step-2 deflation / depth tower

**Date:** 2026-06-10 (code session)
**Engine:** `jobB_tower.py`, `jobB_probe.py`, `jobB_probe2.py`, `core_quotient.py`.
**Recall:** for each tie (|J\*|≥2) the leading `π^μ` coefficient of `G_λ(i)` cancels; the survivor
depth is `d = v_π(G) − μ` (finite except `(2,2)`, where `d=∞`). Job 3 (2026-06-09) found the naive
Step-2 reduction FALSE — 918/1624 ties need the `val>μ` tail. Job B charts that tail.

## Method (carry-honest)

`term_j = C(m,j) iʲ (1+i)ʲ M_j`, with `v_π(term_j)=val(j)`. Process terms by **ascending
val-level** `v₀=μ < v₁ < …`; let `P_t = Σ_{val(j)≤v_t} term_j` (exact ℤ[i]) and record `v_π(P_t)`.
This is the true deflation trace — no mod-2 digit guessing. Also:
- `d₁ = v_π(Σ_{j∈J*} term_j) − μ` (within-J\* / first-level depth),
- the **surviving-order index set** = `{j : digit of u_j at position v_π(G)−val(j) is 1}`.

Output: `results/step2-tower.csv` (1624 ties, m≤12): `λ, m, μ, vpi, d, d1, surv_set, surv_is_box,
core4, quot4, parity, v2f, nlevels`.

## What the tower looks like (the worked cases)

```
λ=(3,3,1,1) m=4   μ=6  J*={0,4}   v_π(G)=8   d=2   d1=4
   level 6:  +{0,4}  -> v_π=10        (within-J* cancels to depth 4)
   level 7:  +{1,3}  -> v_π= 8        (next level LOWERS the survivor to depth 2)
   level 12: +{2}    -> v_π= 8
   surviving-order set {0,1}  (a box, but see below — atypical)

λ=(6,6,4,4,2,2) m=12  μ=20  J*={0,4,8,12}  v_π(G)=42  d=22  d1=14   (largest d)
   level 20: +{0,4,8,12} -> v_π=34     (within-J* depth d1=14)
   level 23: +{1,3,9,11} -> v_π=26
   level 24: +{2,10}     -> v_π=26
   level 25: +{5,7}      -> v_π=30
   level 30: +{6}        -> v_π=42     (TOP level RE-RAISES v_π by 12!)
   surviving-order set {0,4,6,8,10}  (NOT a box)
```

The second case is the signature of the whole tail: **a multi-level cascade**, where even the very
last val-level (`{6}`, val 30) cancels against the accumulated remainder and pushes `v_π` *up* by 12.
Step 2 is not a single deflation — it is a chain of them.

## The three asks, answered

### (1) Is the surviving-order index set a second-level box?  **NO.**

`surv_is_box` true on only **469/1623** ties (≈29%). The largest-d case `(6,6,4,4,2,2)` has surviving
set `{0,4,6,8,10}` — not affine-2-adic. The "deflate to a smaller box with its own unique-min" hope
**fails**.

### (2) Is `d` predicted by the 4-core / 4-quotient?  **NO (only vacuously, via λ).**

| predictor | single-valued |
|---|---|
| 4-core alone | 34/39 — **misleading** |
| 4-quotient alone | 252/416 (≈61%) |
| `\|J*\|` alone | 0/2 |
| `d₁` (within-J\* depth) | 2/13 |
| `(4-core, 4-quotient)=λ` | 1623/1623 — **vacuous** |

The "34/39" is *not* a clean law: the 5 multi-valued cores carry **611/1623 ties**, including the
dominant 4-core `(2,2)` (415 ties, `d∈{2,4,5,6,7,8}`) and three more with n=59. Several
single-valued cores do have real mass (four with n=164), so the 4-core *partially* organises `d` —
but it does **not** determine it. Within a fixed core the 4-quotient finishes the job (e.g. core
`(2,2)`: quotient ⟹ d, 415/415) — but that is just `(core,quotient)=λ` again. **No simple quotient
statistic works**: within core `(2,2)`, `d` is not a function of `|quotient|`, of the number of
nonempty quotient parts, or of `max|q_r|`.

> **Partial pattern (the one clean shard):** within 4-core `(2,2)`, **odd `|quotient|` ⟹ d=2**;
> even `|quotient|` gives the spread `{4,5,6,7,8}`. A parity gate on the quotient size — worth a
> closer look, but far from a full predictor.

### (3) The within-J\* sum is usually NOT the survivor.

`d` vs `d₁` over 1623 ties:
- `d = d₁` (within-J\* sum is the survivor): **705**
- `d < d₁` (next levels LOWER the survivor): **218**
- `d > d₁` (within-J\* over-cancels; next levels recover at higher order): **700**

So **918 = 218+700** ties need the `val>μ` tail — exactly Job 3's count. The 700 `d>d₁` cases are
the cascade: the within-J\* sum cancels *too far*, and later levels cancel against its remainder to
lift `v_π` back up. `d₁` predicts `d` only 2/13 — the first level does not control the survivor.

### Cleanest correlate: `d` grows with `|J*|`

| `\|J*\|` | min d | max d | distribution |
|---|---|---|---|
| 2 | 1 | 8 | `{1:484, 2:702, 4:74, 5:24, 6:53, 7:26, 8:2}` |
| 4 | 4 | 22 | `{4:115, 5:12, 7:64, 9:38, 11:14, 12:14, 22:1}` |

Bigger box ⟹ deeper tower. This is the only robust monotone trend.

## Verdict — what this hands PROVE

The optimistic Step-2 picture (deflate to a *second-level box* with its own unique-min, recurse) is
**refuted**: the surviving-order set is not a box (29%), and `d` is governed by the full λ, not by the
4-core or 4-quotient or any simple quotient statistic. The tail is a genuine **multi-level π-adic
cascade** in which late val-levels actively re-raise `v_π` (700/1623 over-cancel-and-recover).

**Implication:** Step 2 should NOT be attacked as a recursive box reduction. The right tool is a
**whole-series π-adic argument** — e.g. the exact lift `Φ(z)` / generating function (cf. the d=4
engine memory) summed over ALL levels at once — rather than peeling J\* then a second box. The one
foothold worth pursuing: the **parity gate** (core `(2,2)`, odd `|quotient|` ⟹ d=2) and the monotone
`d`-vs-`|J*|` growth, which together hint the depth is controlled by the box dimension `log₂|J*|`
plus a quotient-parity correction — not by a self-similar sub-box.
