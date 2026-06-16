# FINDINGS — Job B (2026-06-16 code session): general-`c` subtop scout

**Question.** Does the `c=3` boundary subtop structure
(`M_{b+c−i} = (poly in a,b)·N_i/const`) persist uniformly in `c`?

**Verdict (nuanced).** The **top** and **first subtop** persist *cleanly and uniformly for
all `c`* — I derived and verified their closed forms for `c=2,3,4,5`. But the **deeper
subtops break the c=3 pattern at `c=4`**: the second subtop `M_{b+c−2}` factors into two
a-linear pieces *only* for `c=3`; for `c≥4` it is an **irreducible quadratic in `a`**. So a
single "multiply in one a-linear factor per descent step + Lemma F2" recipe closes the top
**two** boundary indices for every `c`, but **not** the full boundary for `c≥4`. The
mechanism (direct `Δ` + factor-in-product) is uniform; the *bookkeeping* is not.

All forms fitted exactly from the signed three-row Jacobi–Trudi `M_j` (`mn.py`) and
factored over `ℤ[a,b]`. Script: `threerow-boundary/generalc_subtop_scout.py`.

---

## 1. TOP and FIRST SUBTOP — uniform in `c` (verified `c=2,3,4,5`)

**Top (Lemma T, a-FREE hook):**
> `M_{b+c} = (b−c+1)·∏_{t=2}^{c}(b+t) / c!`

**First subtop (one a-linear deficit):**
> `M_{b+c−1} = (b−c+1)·∏_{t=2}^{c−1}(b+t)·N_{c−1}^{(c)} / c!`,
> with **`N_{c−1}^{(c)} = a(b+c) − (b² + (c−1)b + c(c−2))`**.

Verified (fit == conjectured form, symbolic zero) for `c=2,3,4,5`:

| `c` | `N_{c−1}^{(c)}` |
|---|---|
| 2 | `a(b+2) − (b²+b+0)` |
| 3 | `a(b+3) − (b²+2b+3)` |
| 4 | `a(b+4) − (b²+3b+8)` |
| 5 | `a(b+5) − (b²+4b+15)` |

(constant `c(c−2)` = 0, 3, 8, 15 — exact.) The descent from the top replaces the lone
`(b+c)` factor by `N_{c−1}^{(c)}`, whose leading a-term `a(b+c)` recovers the hook×a.

## 1b. SECOND SUBTOP — where uniformity BREAKS at `c=4`

`M_{b+c−2}` carries two a-factors. Factoring the deficit bracket over `ℤ[a,b]`:

| `c` | second-subtop deficit | a-structure |
|---|---|---|
| 3 | `(a−b+1)(ab²+5ab+6a−b³−b²−4b−6)` | **two a-linear** (factor-in-product works) |
| 4 | `a²b²+7a²b+12a²−2ab³−9ab²−21ab−20a+b⁴+2b³+13b²+16b−8` | **irreducible a-quadratic** |
| 5 | `a²b²+9a²b+20a²−2ab³−13ab²−45ab−70a+b⁴+4b³+29b²+56b+30` | **irreducible a-quadratic** |

So the c=3 "each descent step multiplies in one *linear* `a`-factor `N_i`" picture is
**c=3-special**. For `c≥4` the deeper deficits are irreducible higher-degree-in-`a`
polynomials — they are *not* products of `(P+t)`-type members, so the factor-in-product
engine does not apply term-by-term to them.

## 2. The a-even `v₂(N)` shortcut is also c=3-special

The c=3 proof leans on **`v₂(N₂)=1` for a even**. The first-subtop deficit's a-even
2-adic distribution is **not** uniform:

| `c` | a-even `v₂(N_{c−1}^{(c)})` distribution |
|---|---|
| 2 | spread `{1:777, 2:370, 3:185, …}` |
| 3 | **`{1: 1442}`** (constant — the special case) |
| 4 | spread `{1:684, 2:378, …}` |
| 5 | `{2:990, 3:171, …}` (min **2**, never 1) |
| 6 | spread `{1:665, 2:315, …}` |

Only `c=3` has the constant `v₂=1`. So the clean a-even shortcut does **not** generalize;
a general-`c` proof must keep the deficit attached to its consecutive product (direct `Δ`).

## 3. Boundary lemma still holds for `c=4,5` (margin), and F2 is c-independent

- Direct val sweep (`m ≤ 42`, crude box gate `b≥2c`): `Δ(b+i) > −θ` with **0 violations**;
  min `Δ` = **8** (`c=4`, at `(66,10,40)`), **1** (`c=5`, at `(59,12,38)`). The boundary
  loses for `c=4,5` too — the *phenomenon* is uniform even though the proof tooling is not.
- **Two-factor Lemma F2 is genuinely c-independent**: it is a statement about consecutive
  integers only (`Q≥6`, remove any two designated members, the rest is even). Verified
  **0 failures over 1 616 000 checks** (`Q=6..29`, `R=0..399`, all pairs). It is the same
  keystone for every `c`; the c-dependence is only in *which* members are designated and in
  `Q = b//2 + offset` (which exceeds 6 throughout each box interior).

---

## On-ramp assessment for a general-`c` three-row boundary theorem

**What transfers for free (all `c`):**
- Closed forms for the **top** and **first subtop** (Lemma T + `N_{c−1}^{(c)}`).
- The **direct-`Δ`** method (positive-linear + Kummer carries − 2·deficit).
- The **F2 keystone** (c-independent), closing the top two boundary indices.

**The genuine obstruction (`c≥4`):**
- Deeper boundary indices `b+i`, `i ≤ c−2`, carry deficits that are **irreducible
  polynomials of a-degree `c−i`**, not products of `(P+t)` members. The factor-in-product
  argument needs a *new handle* for these — candidates: (a) a `k`-factor Lemma F_k applied
  after a Newton-polygon / discriminant split of the irreducible deficit, or (b) bounding
  `v₂` of an irreducible a-quadratic via its `a mod 2^k` residue (the a-even/a-odd split
  generalized). Neither is in hand.

**Bottom line.** The c=3 boundary is closed and the **structure** is uniform, but the
**closure technique is not yet uniform past the top two indices**. The on-ramp reaches the
first subtop for all `c`; the interior `NL_c` is already proved
(`proofs/2026-06-14-numberlemma-general-c.md`), so the remaining gap for a general-`c`
theorem is precisely the **deeper boundary deficits** above. Recommend the dream/PROVE
target be: a `k`-factor lemma (or residue bound) for irreducible a-quadratic/cubic deficits.

## Files
`threerow-boundary/generalc_subtop_scout.py` (forms `c=2..5`, a-even `v₂` distributions,
`c=4,5` boundary sweep, F2 c-independence), `boundary_forms.py` (base fitter), `mn.py`.
