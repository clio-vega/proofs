# Two-row d=4 fiber law as an imaginary-part non-vanishing

**Date:** 2026-06-06 (prove session)
**Full write-up:** `2026-06-06-tworow-d4-imaginary-reduction.tex` (compiles, 4pp)

## Theorem (target)
For λ=(a,b), a+b=2m, a≥b≥0: `G_{(a,b)} = ⟨s_{(a,b)}, ψ^m⟩ = 0 ⟺ (a,b)=(2,2)`, ψ=h₂+ie₂.

## What is proved here (the reduction)
Writing b for the second row, `G_{(2m-b,b)} = [u^b]((1−u)(1+su+u²)^m)`, s=1+i, and
```
Im G_{(2m-b,b)} = Σ_{k odd ≥1} (−1)^{(k−1)/2} C(m,k) [ T(m−k, b−k) − T(m−k, b−k−1) ]
```
with `T(N,j)=[u^j](1+u+u²)^N` the trinomial coefficients (≥0, unimodal). Each bracket is ≥0
(trinomial unimodality), so **Im G is an alternating sum of nonnegative integers** c₁−c₃+c₅−…,
last index ≤ b.

**Reduction (Prop. 1):** Im G ≠ 0 ⟹ G ≠ 0. So the law follows from
> **Im G_{(2m-b,b)} ≠ 0 for all 1≤b≤m, (m,b)≠(2,2).**

This is a *discrete, robust* statement — unlike the earlier continued-fraction
reduction whose margin looked like a knife-edge 1/4 (that was a normalization artifact of
dividing by n).

## Proof — complete for boundary families b=1,2,3,4
- **b=1 (hook):** `G_{(2m-1,1)} = (m−1) + m·i` (direct: r₀=1, r₁=ms). Im G = m. |G|²=m²+(m−1)².
- **b=2:** `G_{(2m-2,2)} = m(m−2)·i` (Re G = 0). **The unique vanisher**, at m=2 only = (2,2).
- **b=3:** Im G = m(m−1)(m−2)/3 ≠ 0 for m≥3.
- **b=4:** Im G = m(m−1)(2m−7)/3 ≠ 0 for integer m≥4 (2m−7 odd).

## Sharp structure (verified, not yet proved)
- **min over b of |Im G| = m, attained uniquely at the hook b=1** (m≥4); = 2 at m=3 (b=3).
  So `|Im G_{(2m-b,b)}| ≥ m` — the law is far from marginal.
- Verified exhaustively: Im G ≠ 0 except (2,2), and |Im G| ≥ m, for **all 3≤m≤300, 1≤b≤m**.
- Integer roots of I_b(m):=Im G (poly in m) are exactly {0,1,…,⌊(b−1)/2⌋}, all < b, hence
  outside the valid range m≥b. Every root of I_b in [b,∞) is irrational.

## Gap (precise)
Uniform non-vanishing of the alternating trinomial sum for **b≥5** — equivalently, for each b,
the "tail" factor of I_b(m) (after the forced linear factors m,m−1,…,m−⌊(b−1)/2⌋) has no
integer root ≥ b. Per-b this is Diophantine; the difficulty is uniformity in b.

Two concrete routes (see §5 of tex):
- **Route A (2-adic):** C(m,k) and T(N,j) mod 2 are digit-governed (Lucas + trinomial
  self-similarity); the leading 2-adic term of the sum should be computably nonzero
  ⟹ finite v₂(Im G). [NB: Im G is even for all even m, so this is subtle.]
- **Route B (ballot):** Im G is a signed count over length-m words in {0, real-1, imag-1, 2}
  weighted i^{#imag}; a sign-reversing involution with nonempty fixed-point set gives ≠0,
  and is the natural candidate to also yield the sharp |Im G| ≥ m.

## Verification scripts
`~/projects/scratch/2026-06-06-prove/`: rl.py, closed.py, imstruct.py, improfile2.py, verify_big.py.
