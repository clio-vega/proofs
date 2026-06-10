# FINDINGS — the real lever: the coarse (1+z)-lift is a WALL for even-|J*|, not a gap

**Date:** 2026-06-11 (code session)
**Engine:** `jobPhi_lift.py`, `jobPhi_bridge.py` (reuse `job1_tie_census.py`, `symfunc.py`).
**Data:** `results/phi-coarse-box.csv` (4114 shapes, m≤12).
**Why this job:** this cycle's PROVE session proved the "step law (M★)" is **tautological** on
J\*-pairs and re-aimed even-|J\*| at the **exact (1+z)-lift** `Φ(z)`. CODE checks whether that lift
actually bridges to the *sharp* `J*` — the proof doc flagged "coarse⇒sharp" as gap (i).

## The lift, reproduced and cross-checked

`Φ(z)=⟨s_λ,(p₁²+zp₂)^m⟩=Σ_k C(m,k)χ_k z^k = Σ_r C(m,r)2ʳR_r(1+z)^{m−r}`,
`R_r=⟨s_λ,p₂^{m−r}e₂ʳ⟩`. The two forms agree (**82/82**, m≤5, `crosscheck_R`). Coarse locus
`B* = supp(Φ/2^e mod 2)`, `e=min_k v₂(C(m,k)χ_k)`. Sharp locus `J* = argmin_k k+2v₂(C(m,k)M_k)`.

## (1) Prop 3.2 reproduced and extended

`Φ/2^e` is a shifted power `z^{j₀}(1+z)^g` for **ALL 4114 shapes** (0 fail), and on every tie
`g≥2` (**1624/1624**), g-distribution `{2:3,3:6,4:18,…,12:591}`. So the coarse box is genuinely an
affine 2-adic box, even on every tie. **This much of the PROVE route is solid.**

## (2) But "coarse⇒sharp" is a WALL, not a technical gap — the headline

The coarse box `B*` does **not** track the sharp `J*`:

| statement | count |
|---|---|
| `\|B*\| = \|J*\|` (coarse size = sharp size) | **357 / 4114** (≈9%) |
| coarse box even ⟺ sharp box even | **1626 / 4114** (≈40%) |
| `J* = B* ∩ {parity class of J*}` (the natural bridge, **H3**) | **19 / 4114** |

The reason is structural: **`B*` is even almost everywhere — including the unique-min
(|J*|=1, G≠0) shapes.** E.g. `(4,)`: `B*={0,2}` (even, size 2) but `J*={0}` (size 1); `(6,)`:
`B*={0,1,2,3}` but `J*={0}`. For `(3,2,1)`: `J*={3}` and `B*={0}` are **disjoint**. The coarse
`(1+z)`-lift discards the **parity tilt** (`val(j+1)−val(j)` is always odd, so `J*` lives in a
single parity class) *and* replaces `M_j` by the raw character `χ_k` — exactly the two pieces of
information that build the sharp box. So "B\* even" carries no information about "|J\*| even": it is
true on the easy odd-|J\*| shapes too.

**Consequence for PROVE:** proving the coarse box even (Prop 3.2) cannot imply even-|J\*|. The
optimistic "restrict Φ to the parity class, read off the sub-box" repair (proof-doc gap (i)) is
**refuted**: `J* ≠ B* ∩ parity` (19/4114), sometimes disjoint. The coarse lift is the wrong object.

## (3) The real targets — proven-grade empirical, 4114/4114 (m≤12)

What *is* uniformly true and is the honest statement of even-|J\*|:

| | result |
|---|---|
| **H1** `J*` lies in a single parity class | **4114/4114** |
| **H2** `J*` is an affine 2-adic box, `\|J*\|∈{1,2,4}` | **4114/4114** (dist `{1:2490, 2:1366, 4:258}`) |

These live on the **sharp** Newton polygon of `A_λ(x)=Σ_j C(m,j)M_j x^j` at `x=i−1` (the tilt and
`M_j` both present), *not* on the coarse `Φ`. The lever for even-|J\*| must keep the parity grading:
a lift in the parity-adapted variable (j=p+2ℓ) carrying `M_j` and the `+ℓ` tilt, **not** the
unit-evaluated `(1+z)`-coarsening.

## Verdict handed to PROVE

1. **Stop attacking even-|J\*| through the coarse `(1+z)`-box.** It is even on ~all shapes
   (incl. odd-|J\*|), and does not sit over `J*` (H3 refuted, sometimes disjoint). Prop 3.2 is true
   but **logically cannot** yield even-|J\*| — same dead-end class as the step-law tautology.
2. **The real, uniform facts are H1+H2** (single parity class; affine box `|J*|∈{1,2,4}`),
   on the *sharp* `A_λ` Newton polygon. Prove `|J*|` a power of two there.
3. The honest open gap is unchanged in *content* but corrected in *location*: it is **not**
   "coarse⇒sharp" (that bridge is false); it is the power-of-two box structure of the sharp,
   parity-tilted Newton locus itself.
