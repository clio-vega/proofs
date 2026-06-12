# Hooks done: `|J*|` even proved for the whole hook family (d=4)

Robin — short version: the 06‑12 prove session closed the `|J*|`‑even question **for every hook**
`λ = (a,1^b)`, completely and rigorously (no gaps). The general‑λ case is still open, but I now know
exactly why, and the hook proof gives a concrete template for the next family.

**The theorem (hooks).** For every hook `λ ⊢ 2m`, `|J*| ∈ {1,2}` — so `|J*|` is even whenever it's
`≥ 2`, and the leading‑π layer cancels on every hook tie. The only possible tie is `J* = {0,2}`
(happens iff `m` odd and `b ≡ 2,3 mod 4`).

**Why it worked (and what's new):**
1. **Hook closed form** `M_j = C(2m−1−j, a−1)` — proved via reproducing kernel
   (`h₂+ce₂` becomes `(s+t)(s+ct)` under `p₁=s+t, p₂=s²−t²`) + dual Jacobi–Trudi + a telescoping
   recursion. Bonus: `G_λ(i) = [z^{a−1}](1+z)^{m−1}(z+i)^m` for hooks (a clean non‑vanishing handle).
2. **Falling‑factorial cancellation** collapses `val(j)−val(0)` to a pure Kummer expression
   `j + 2[v₂C(m,j)+v₂C(b,j)−v₂C(2m−1,j)]`.
3. **Two carry lemmas** (Lemma K / K′) bound `v₂C(m−1,⌊j/2⌋) − v₂C(m,j) ≤ ⌊j/2⌋`, with equality only
   at `j=2`. The proof is two lines of `v₂`‑identities + a parity split — genuinely clean.

**Two things I want to flag (the methodology lessons):**
- The PROVE.md warning was right: I checked that **Newton‑polygon/ramification formalism alone cannot
  prove `|J*|` even** — it only gives "edge length even", which is the trivial single‑parity fact H1.
  Evenness is a real arithmetic property of the `M_j`. So the hook win is *not* a formal trick; it's
  the first place the arithmetic is fully exposed.
- I also nailed down that the **binomial skeleton is trivial** (`argmin (j+2v₂C(m,j)) = {0}` always):
  every tie is created entirely by `v₂(M_j)`. That's the precise shape of the 4‑core obstruction.

**Sharper conjecture (verified `m ≤ 10`, 0 exceptions):** in general the on‑line set is always a
**2‑adic box**, so `|J*|` is a **power of 2** (`∈ {1,2,4}` observed). Even is the weak form; power‑of‑2
is the real structure. I did *not* prove this in general — and per the argmin lesson I deliberately did
not manufacture an optimiser‑internal "law".

**Next.** Two‑row shapes `(a,b)`. They should have a Prop‑2 analogue (trinomial generating function),
which is the natural way to extend the falling‑factorial cancellation. If that works, hooks + two‑rows
covers a lot of the boundary and might suggest the general comparison generating function.

Full writeup + proofs: `projects/proofs/2026-06-12-hook-Jstar-even.md`. Verification script:
`projects/code/2026-06-12-hook-Jstar-even-verify.py` (all checks pass; hooks `m<300` for the main
inequality, Lemma K to `m<2000`).

— Clio
