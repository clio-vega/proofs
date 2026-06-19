# FINDINGS — Job B (2026-06-19 CODE): route-1 q-determinant probe for the general-`c` Content Lemma

**Scripts:** `threerow-boundary/jobB_route1_qdet.py` (q-hook determinant + tower),
`jobB_route1_straighten.py` (bialternant straightening test).
**Builds on:** `FINDINGS-2026-06-18-jobB-content-decomp.md` (§Route-1 stub), 07-05 crown
`2026-07-05-content-lemma-four-routes-converge` (route 1 = highest-leverage).
**Target:** does the boundary deficit `N_i^{(c)}` have a q-analogue whose **full** cyclotomic-tower
multiplicity `Σ_r mult_{Φ_{2^r}}` equals its exact 2-adic content `g[c][k]` (including the `r≥2`
surplus the `q=−1`-only stub missed)?

## Verdict (one line): **NO / PARTIAL — the q-hook (dimension) construction does not reach `M_{b+i}`.**

The **cyclotomic-tower *principle* is confirmed** (a leak-free q-lift's tower = `v₂` exactly), but the
three-row alternant `M_{b+i}` is **not a dimension / q-hook product** — not even in the interior — so
the dimension q-determinant is the wrong q-lift. Route 1 still *could* close, but only with the **LGV
Gaussian-binomial determinant of the closed form**, not the hook formula. This sharpens the browse
target for `2502.06032` and **corrects a latent false belief** (`M_j = f^{(a−j,b−j,c−j)}`).

## What was established

### 1. The cyclotomic-tower principle holds (the valid foundation) — `jobB_route1_qdet.py §0`
For any q-polynomial `P ∈ ℤ[q]`, `Σ_{r≥1} mult_{Φ_{2^r}}(P) = v₂(P(1)) − leak(P)`, where `leak ≥ 0` is
the 2-content hidden in `P`'s non-cyclotomic factors / leading coeff (since `Φ_{2^r}(1)=2`,
`Φ_d(1)=1` for non-prime-power `d`, `=p` odd for odd-prime powers). A q-lift that is a genuine
**ratio of q-integers** (q-hook formula `[n]_q!/∏[h]_q`) has **`leak = 0`**, so its tower multiplicity
equals `v₂` of the `q=1` value exactly. Verified on 7 genuine dimensions
(`(3,2,1),(4,2,2),…,(6,5,4)`): `tower == v₂`, **0 anomalies**. So *if* `M_{b+i}` were a dimension,
route 1 would close trivially.

### 2. But `M_j` is NOT a dimension for three rows (the premise fails) — verified vs MN
The route-1 plan implicitly assumed `M_j = f^{(a−j,b−j,c−j)}` (true and clean for **two** rows:
`M_j = f^{(a−j,b−j)}`). For **three** rows this is **false except at `j=0`**:

| `λ=(10,8,6)`, `m=12` | `j=0` | `j=1` | `j=2` | `j=3` | `j=4` | `j=5` |
|---|---|---|---|---|---|---|
| `M_j` (MN) | 267711444 | 100876776 | 38135318 | 14464008 | 5503914 | 2101125 |
| `f^{(a−j,b−j,c−j)}` | 267711444 | 14284998 | 787644 | 45045 | 2673 | 162 |
| equal? | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ |

(Same pattern for `(3,2,1)`, `(5,4,3)`.) The inner product `M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩` is a
**signed combination of characters**, not a single hook product, once there are ≥3 rows. *(Note: the
header comment of `threerow-c3/mn.py` calls `f^{(a−j,b−j,c−j)}` the "conjecture"; its `test_threerow`
actually checks `M_j` against the abacus-χ engine, not against `fdim`, so nothing downstream relied on
the false identity — but the comment is misleading and should be flagged.)*

### 3. At the boundary the dimension q-determinant literally vanishes (24/24) — `jobB_route1_qdet.py §1`
For the deficit indices `j = b+i` (`i ≥ 1`), `mu = (a−j, b−j, c−j)` has **two negative parts**
(`b−j = −i`, `c−j < 0`). The dimension-determinant continuation `[n]_q! · det(1/[mu_p−p+r]_q!)` is
**`0`** at all 24 tested `(c,k,parity)` cells (`c=4..9`, `k=4,5,6`), while `M_{b+i} ≠ 0` (e.g.
`c=5,k=4`: `D_q=0`, `M=332406`). The bialternant **does not straighten** to a partition
(`jobB_route1_straighten.py`: `0/24` give a valid `ν`) — confirming `M_{b+i}` is an alternant
continuation that is *not* `±` any dimension.

## Why this is a useful negative (and what route 1 actually needs)
- The content of `M_{b+i}` is genuinely **cancellation-born** (cf. the stub's `c=7,k=5`: the
  alternating sum lifts `v₂(M)` by `+3` over the minimal term). A construction whose value is a single
  hook product (min-of-parts content) **cannot** see it.
- The correct q-lift is the one the note named: the **LGV / Gatzweiler–Krattenthaler q-determinant of
  the closed form** — the `3×3` Gaussian-binomial determinant from the `D_j` expansion
  (`M_j = Σ_{σ∈S₃} ± ∏ C(N,·)` with each `C(N,k) → [N,k]_q`). By LGV this is **q-positive** (lattice
  paths in a box ⟹ no alternating cancellation), and being a sum of products of Gaussian binomials it
  is a controlled object whose `Φ_{2^r}` tower can be read off. The dimension hook is *not* this object
  (it collapses to 0 at the boundary).
- **What to test next (post-browse):** lift the closed-form `C(N,b−j)`-determinant (not the hook) and
  check `Σ_r mult_{Φ_{2^r}}` vs `g[c][k]` on the same grid. The tower *principle* (item 1) guarantees
  that IF that q-determinant is a ratio of q-integers (or LGV-positive with no 2-leak), the match is
  exact; the open question is whether the LGV determinant has `leak = 0` or whether the
  cancellation-born surplus survives as non-cyclotomic 2-content.

## Honesty / scope
- The cross-check `D_q(1)` vs `M` (the q→1 gate) was run on **every** case and **failed everywhere at
  the boundary** (D vanishes) — reported loudly rather than patched. This killed the dimension-hook
  route cleanly.
- A clean NEGATIVE with a structural reason (cf. Ayyer–Kumari): route 1 is **not closed** by the
  dimension q-hook; the remaining gap is the **q-determinant construction** (LGV Gaussian-binomial of
  the closed form), not a cyclotomic-tower mismatch. The tower principle itself is sound and verified.
