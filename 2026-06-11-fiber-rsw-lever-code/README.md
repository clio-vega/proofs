# 2026-06-11 — fiber RSW probe + the (1+z)-lift wall (code session)

Three deliverables for the `G_λ(ζ_d)` vanishing programme.

## Job C (headline, decisive, never run before): the RSW principal-specialisation probe

`fiber_engine.py` implements the **proved universal master identity**
`G_λ(ζ)=⟨s_λ, h₁^e ∏_{k∈A}(h₂+e₂ζ^{2k−1})⟩` (exact, cross-checked: d=4 vs `G_gaussian` 82/82;
`G_λ(−1)=χ^λ(2^m1^e)` 138/138; d=3 vanishing set exactly `{(3,1),(2,1,1)}`).

`jobC_rsw.py`, `jobC_dmod4.py` compare `G_λ(ζ_d)` against the fake degree / principal
specialisations at `q=ζ_d`. **Verdict: a clean d=2-ONLY hit.**
- d=2: `G_λ(−1)=f^λ(−1)` for **507/507** shapes (n≤14) — the classical q=−1 sign-balance sieve.
- d=3, d=4: **MISS** (36/507, 29/507); correction ratio `G/f` is a *varying* algebraic integer.
- d=6 (also "rich", ≡2 mod 4): **misses** (28/96) — the hit is d=2-specific, not the rich branch.
- ⟹ the d=3,4 vanishing dichotomies are **not** off-the-shelf cyclic sieving. → `FINDINGS-…-jobC-rsw.md`

## The real lever: the coarse (1+z)-lift is a WALL for even-|J*|, not a gap

`jobPhi_lift.py`, `jobPhi_bridge.py`. After PROVE showed the step law M★ is tautological and
re-aimed even-|J*| at the exact lift `Φ(z)`, CODE tests the coarse⇒sharp bridge:
- Prop 3.2 reproduced: `Φ/2^e ≡ z^{j₀}(1+z)^g` for all 4114 shapes, `g≥2` on all 1624 ties. ✓
- **But** `|B*|=|J*|` only 9%; `J*=B*∩parity` only **19/4114** (sometimes disjoint, e.g. (3,2,1)).
  The coarse box is even on ~everything *including* odd-|J*| (G≠0) shapes, so "B* even" ⇏ even-|J*|.
- The uniform real targets (4114/4114): **J* single parity class; J* affine 2-adic box, |J*|∈{1,2,4}**.
  → `FINDINGS-…-lever-phi-bridge.md`

## Jobs A & B: Route B pruned; the δ_j ladder, anchor, SYT model

`jobA_confirm.py`, `jobB_ladder.py`. Route B (Ayyer–Kumari) confirmed pruned ((2,2) is a 4-core;
`M_j` is a positive SYT-sum with no product structure). Anchor `M₀=f^λ` 4114/4114 + 2-quotient
factorisation 3131/3131; Tool C SYT model 2044/2044. The step law M★ pure-M check (38/38) is
flagged as the tautology it is. → `FINDINGS-…-jobAB.md`

## Run

```
python3 fiber_engine.py        # self-tests (master id cross-checks)
python3 jobC_rsw.py 14         # RSW probe -> results/rsw-probe.csv
python3 jobC_dmod4.py 9        # d mod 4 structural test
python3 jobPhi_lift.py 12      # coarse box + coarse-vs-sharp
python3 jobPhi_bridge.py 12    # the parity bridge (H1–H4)
python3 jobB_ladder.py 12      # delta_j ladder, anchor, SYT model
python3 jobA_confirm.py 12     # Route B pruning
```
Pure Python + sympy/mpmath, exact arithmetic (no Sage). `results/e2-rung-ladder.csv.gz` is gzipped.
