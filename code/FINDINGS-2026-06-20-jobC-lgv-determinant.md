# FINDINGS — Job C (2026-06-20 CODE): the LGV Gaussian-binomial q-determinant of `M_j`

**Script:** `jobC_lgv_qdeterminant.py` · **Ground truth:** MN (`threerow-c3/mn.py`),
determinantal form from `threerow-c1/dj.py`.
**Builds on:** the 06-19 route-1 negative (`FINDINGS-2026-06-19-jobB-route1-qdet.md`) — the dimension
q-hook **vanishes** at the boundary; the correct q-lift is the **Gaussian-binomial determinant of the
closed form**, which this job implements and tests.

## The object

`M_j = Σ_{σ∈S₃} sgn(σ) D_j(L₀+σ₁, L₁+σ₂, L₂+σ₃)`, `L=(a−1,b−2,c−3)`,
`D_j(p,q,r) = [s^p t^q u^r] e₁^{2(m−j)} e₂^j = Σ_{i+k+l=j} multinom(j;i,k,l)·multinom(N; p−i−k, q−i−l, r−k−l)`,
`N=2(m−j)`. Each multinomial is an LGV box-path count, so `M_j` is a signed `S₃` sum of products of
ordinary binomials. The **q-lift** replaces every `multinom(·)` by its **Gaussian q-multinomial**
(nested Gaussian binomials = inversion-statistic box-path generating function).

## Results — the three plan questions

### (1) q=1 recovery — and it REACHES THE BOUNDARY  ✅ (the key advance)
`q=1` value of the q-determinant **= `M_j` exactly: 35/35 cases, 0 mismatches**, over interior *and*
boundary indices `j = b+1, b+2`. This is exactly where the **dimension q-hook collapsed to 0** (06-19,
24/24). So the Gaussian-binomial closed-form determinant is the **right q-lift**: it is a genuine
nonzero q-object at the deficit indices `M_{b+i}`.

### (2) q-positivity  ✅
The q-lift is **q-positive coefficient-wise in every case: 27/27, 0 negative coefficients** — despite
being a signed `S₃` sum. This confirms the LGV prediction: the alternating cancellation telescopes to a
nonnegative box-path count (the inversion-weighted q-multinomial carries the right statistic). The
naive Gaussian q-multinomial lift already suffices for positivity here; no ad-hoc q-power correction
was needed.

### (3) Leak-free?  ❌ NO — the 2-content escapes the cyclotomic tower
Reading `v₂(M_{b+i})` off the `Φ_{2^r}` (=`q^{2^{r-1}}+1`) tower of the q-positive determinant
**FAILS**: the cancellation-born 2-content survives as **non-cyclotomic** content.

```
 λ          j=b+i   M_j     v₂(M)  Φ_{2^r}-tower  leak (= v₂ − tower)
 (8,6,4)    7       1059    0      0              0
 (8,6,4)    8       414     1      1              0      <- clean (content in q^?+1)
 (10,8,6)   9       45936   4      0              4      <- LEAK: all 4 outside the tower
 (10,8,6)   10      17769   0      0              0
```

At `(10,8,6) j=9` the full `v₂ = 4` lives in **non-2-power-cyclotomic factors** (a q-positive
polynomial can have even `P(1)` with no `q+1` or `q^2+1` factor). So **`Σ_r mult_{Φ_{2^r}} ≠ v₂` in
general** for this q-lift.

## Verdict for route 1 (the general-`c` Content Lemma)

- **Positive half confirmed:** the LGV Gaussian-binomial determinant of the closed form is the correct,
  boundary-reaching, **q-positive** lift of `M_{b+i}` — the object the 06-19 negative pointed to. It
  exists and behaves as LGV predicts.
- **But ONE positivity theorem does NOT close the Content Lemma** via cyclotomic-tower counting: the
  2-adic content is **not leak-free**. The `v₂` we need is partly (sometimes wholly) carried by the
  polynomial's **non-cyclotomic** factors, which a `Φ_{2^r}`-tower reading cannot see. Route 1 therefore
  needs more than "LGV-positive ⟹ tower = v₂": it needs an account of the non-cyclotomic 2-content
  (e.g. a charge/normalisation that pushes all 2-content into `q^{2^{r-1}}+1` factors — the precise
  Gatzweiler–Krattenthaler quotient class, **arXiv 2502.06032**, still to be pulled).

## Honesty / scope
- q=1 gate run on **every** case (incl. boundary) — 0 mismatches, reported as the trust anchor.
- The leak is demonstrated on a concrete witness `(10,8,6) j=9` (`v₂=4`, tower=`0`), not asserted.
- Caveat: this is the **natural inversion-statistic** q-multinomial lift. A different q-power
  normalisation could redistribute content between cyclotomic and non-cyclotomic factors; the canonical
  G–K normalisation (2502.06032) may or may not be leak-free — that is the open question this probe
  sharpens, not settles. Coefficient `Φ_{2^r}` detection covers `q+1, q²+1, q⁴+1, …` (all 2-power
  cyclotomics); the leak is genuine, not a missed factor.
