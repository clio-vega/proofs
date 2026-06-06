# d=4 imaginary-part reduction — does it reach Gap A? (code session 2026-06-06)

**Verdict: No.** The imaginary-part reduction `G_λ(i)=0 ⟹ Im G ≠ 0` is genuinely
two-row-special. For general `λ ⊢ 2m`, `Im G_λ` vanishes on a *rich* set (70 shapes for
m ≤ 12, of which 52 have `Re ≠ 0`), so it cannot isolate `(2,2)`. The route to Gap A
remains the joint `Re = Im = 0` analysis / 4-core valuation decomposition.

See **`FINDINGS-imaginary-reduction.md`** for the full writeup.

## Files

- `dfour_imaginary_reduction.py` — verifies the all-λ identity
  `G_λ(i) = Σ_k C(m,k) i^k N_{λ,k}` with `N_{λ,k}=⟨s_λ, h₂^{m−k}e₂^k⟩ ≥ 0`
  (three-way cross-check: direct ⟨s_λ,ψ^m⟩ / chain model / bracket sum), the two-row
  trinomial specialization, the `Im G_λ` zero-set (m ≤ 12), and the `|Im G_λ|` bound.
- `dfour_imaginary_structure.py` — conjugation law `G_{λ'}=i^m·conj(G_λ)`, the
  trivial/nontrivial split of the Im-zero set, the full-vanishing anchor
  (`G_λ=0 ⟺ λ=(2,2)`, m ≤ 12), and the two-row sharp bound `min_b |Im G_{(2m−b,b)}|=m`.
- `core.py` — Pieri-rule symmetric-function machinery (self-contained copy).
- `results/` — archived console output of both scripts.

## Reproduce

```
python3 dfour_imaginary_reduction.py     # ~minutes (chain model to m=6)
python3 dfour_imaginary_structure.py     # fast
```

Pure Python + stdlib (no SageMath).
