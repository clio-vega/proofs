# FINDINGS — Job 3: the surviving higher-order coefficient (Step 2)

**Date:** 2026-06-09 (code session)
**Script:** `job3_survivor.py`  **Data:** `results/job3-survivor.csv`

## Step 2 conclusion holds: every tie ≠ (2,2) has finite v_π(G)

After the leading `π^μ` coefficient cancels (`|J*|` even), the surviving cancellation
order `d := v_π(G) − μ` was computed exactly for all 1624 ties (m≤12):

```
d : 1→484  2→702  4→189  5→36  6→53  7→90  8→2  9→38  11→14  12→14  22→1   inf→1
```

- The **lone `inf`** is `(2,2)` — confirming yet again it is the unique total-vanisher.
- `d` is **not constant**: it ranges 1..22. So the deflated value `G_λ(i)/π^μ` carries a
  shape-dependent extra valuation; non-vanishing is robust but the *depth* is irregular.

## Where the survivor comes from — NOT purely within J\*

Let `L₀ = Σ_{j∈J*} u_j(−i)^{v_j}iʲ` be the J\*-only leading unit sum (its v_π measures
how far the within-J\* cancellation reaches). Classifying:

- **705/1624** ties: `d == v_π(L₀)` — the survivor is **within J\***.
- **918/1624** ties: `d > v_π(L₀)` (or L₀ itself vanishes) — the survivor **needs the
  next-nearest-val indices** outside J\*.

## Verdict on the reduction conjecture

**The naive inductive engine — "the deflated `G/π^μ` is governed by a new Newton
polygon on J\* with a unique min" — is FALSE.** In the majority (57%) of ties the
surviving term is not produced inside J\* at all; the `val(j) > μ` tail contributes at
the surviving order. Step 2 therefore cannot be closed by iterating the Newton polygon
on the minimiser locus alone — a correct Step 2 must track the full π-adic tail (or find
a different invariant than "re-deflate and re-minimise"). This is the honest obstruction
for PROVE: the even-|J\*| cancellation (Job 2) is clean and structural, but the *amount*
of surviving valuation is governed by the whole expansion, not by J\*.
