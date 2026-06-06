# Does the imaginary reduction reach Gap A? — VERDICT: **No.**

*Code session, 2026-06-06. Scripts: `dfour_imaginary_reduction.py`, `dfour_imaginary_structure.py`.
Outputs archived in `results/`.*

## The question (IR-4)

The two-row d=4 law had a clean reduction: `G_{(2m−b,b)}(i)=0 ⟹ Im G ≠ 0`, with
`Im G` a signed sum of nonnegative trinomial coefficients. **Does this imaginary-part
reduction generalize to all λ ⊢ 2m, and thereby close Gap A (the full d=4 fiber law
`G_λ(i)=0 ⟺ λ=(2,2)`)?**

Here `G_λ(i) = ⟨s_λ, ψ^m⟩`, `ψ = h₂ + i·e₂ = s_2 + i·s_{11}`, `λ ⊢ 2m`.

## What is true (the good structure)

**1. The signed-nonnegative identity holds for EVERY shape.** Because `s_2` and `s_{11}`
commute,
```
G_λ(i) = Σ_{k=0}^m C(m,k) i^k N_{λ,k},   N_{λ,k} := ⟨s_λ, h₂^{m−k} e₂^k⟩ ≥ 0,
```
so
```
Re G_λ = Σ_{k even} (−1)^{k/2}     C(m,k) N_{λ,k}
Im G_λ = Σ_{k odd}  (−1)^{(k−1)/2} C(m,k) N_{λ,k}            (★)
```
**Verified three independent ways** (direct ⟨s_λ,ψ^m⟩ symbolic expansion; the chain
model; the bracket sum (★)) for **all** λ ⊢ 2m, m ≤ 11 (chain model to m ≤ 6, the
direct+bracket cross-check to m ≤ 11). `N_{λ,k} ≥ 0` throughout. The chain model's
v-distribution satisfies `vdist[k] = C(m,k)·N_{λ,k}` (every interleaving of the same
multiset of strips contributes the same coefficient — commutativity made combinatorial).

The two-row case is the rank-2 instance: `N_{(2m−b,b),k} = T(m−k,b−k) − T(m−k,b−k−1)`
where `T(p,j)=[x^j](1+x+x²)^p` is the central trinomial coefficient. **Verified m ≤ 10.**

So (★) *is* the natural Gap-A generalization of the trinomial formula, and it is a genuine
signed sum of nonnegatives for all λ. `N_{λ,k}` is an iterated-Pieri / domino-type count
(number of chains ∅ ⊂ ⋯ ⊂ λ using m−k horizontal and k vertical 2-strips).

## What is false (the verdict)

**2. `Im G_λ = 0` does NOT single out (2,2).** The imaginary part vanishes on a *rich*
set — **70 shapes for m ≤ 12**, not one:

- **Trivial zeros:** every single row `(2m)` and single column `(1^{2m})` has `G_λ` real
  (`= ±1`), so `Im = 0`. (For a single row, only the all-horizontal `k=0` term survives.)
- **52 nontrivial zeros** (`Re ≠ 0`, `Im = 0`) for m ≤ 12. Examples:
  `(2,2)` [m=2], `(2,2,1,1)` [m=3], `(4,2,2),(4,2,1,1),(3,3,2),(3,3,1,1)` [m=4],
  `(4,4,4,4)` [m=8], `(7,7,2,2,2,2,2)` [m=12], …
  A visible infinite sub-family: `(2m−4j, 2, 1^{4j−2})` and its conjugates
  (e.g. `(10,2,1,1),(6,2,1⁶),(2,2,1¹⁰)` at m=7).

Consequently `min_{λ≠(2,2)} |Im G_λ| = 0` for every m — the imaginary part alone can
never isolate (2,2). **The imaginary reduction is genuinely two-row-special**, exactly as
the earlier Galois-norm proxy was.

## The anchor that still holds

**3. Full vanishing `G_λ(i) = 0` (i.e. `Re = Im = 0` simultaneously) occurs ONLY at
`(2,2)` for m ≤ 12.** Confirmed. So the d=4 fiber-law conjecture stands — but on the
52+ shapes where `Im = 0`, it is `Re G_λ ≠ 0` that does the work. The bridge to Gap A
therefore remains the **4-core valuation decomposition**, not the imaginary part.

## Supporting facts proved/verified

- **Conjugation law:** `G_{λ'} = i^m · conj(G_λ)` — verified all λ, m ≤ 11. Hence the
  *full* vanishing set is conjugation-closed (`(2,2)` is self-conjugate); the *Im*-zero
  set is closed only up to the `i^m` twist (closed for m even, swaps Re↔Im for m odd),
  which is why nontrivial Im-zeros come in conjugate pairs only when m is even (self-conj
  counts: 1,0,2,0,0,0,5,0,0,0,11 for m=2..12).
- **Two-row sharp bound (b ≥ 1):** `min_b |Im G_{(2m−b,b)}| = m`, attained at the hook
  `(2m−1,1)`, for all m ≥ 4. (Exceptions: m=2 → 0 at `(2,2)`; m=3 → 2 at `(3,3)`.)
  This matches the prior "|Im G| ≥ m, minimised at the hook" record and supports PROVE
  Route B for the two-row case.

## Bottom line for the program

- **Keep** the identity (★): a clean, all-λ signed-nonnegative decomposition with
  `N_{λ,k} = ⟨s_λ, h₂^{m−k}e₂^k⟩` a Pieri/domino count. It is the right object.
- **Abandon** the hope that `Im G_λ ≠ 0` alone closes Gap A. It is two-row-only.
- **The route to Gap A is the joint `Re = Im = 0` analysis**, i.e. the 4-core valuation
  decomposition — now confirmed cleanly to be the necessary bridge.
