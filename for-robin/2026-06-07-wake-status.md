# Status for Robin — 2026-06-07 (wake)

Hi Robin,

Daily status. **Gmail is still locked** — the `claude.ai Gmail` connector needs you to run `/mcp`
and authorise it in your browser (this is the ~26th session it's been down). Until then these notes
in `for-robin/` are my only outbound channel. I tried to start the OAuth flow from my side; it just
returns "ask the user to run /mcp and select claude.ai Gmail."

## What last cycle produced

- **PROVE (two-row d=4):** sharpened the last gap to one crisp statement — **(♦) "Q_b has no
  rational root for b ≥ 5"** — and proved real structure around it: a clean reduction `Im G =
  [u^{b−1}]((1−u)H_m)` with `H_m ∈ ℤ[u]` (the `i` fully removed), a forced-root lemma, and the crown
  fact (verified b ≤ 24): `Q_b` is irreducible over ℚ for `b ≢ 0 mod 4`, and `= linear × irreducible`
  when `4 | b`. So the entire two-row law holds for b ≤ 24. Honest negatives: the 2-adic route is
  **dead** (the 2-valuation of `Im G` is unbounded — hits 11 at m=18 — so no finite truncation can
  finish), and the combinatorial involution stalled on a boundary class.
- **CODE:** confirmed the imaginary-part route is **two-row-special** — `Im G_λ = 0` on 70 shapes
  for m ≤ 12, not just (2,2). So the imaginary reduction is NOT the path to the general-λ conjecture.
  Banked the structural identity `G_λ(i) = Σ_k C(m,k) i^k N_{λ,k}` with `N_{λ,k} ≥ 0` (all λ, m≤11).
- **REVIEW:** reviewed Rick's Azenhas–BDI bridge (his P_PARK #5) — endorsed his NEGATIVE verdict,
  but argued he's leaning on the wrong leg: the **3-dim BDI vs 4-dim Azenhas cone** dimension gap at
  n=2 is the airtight obstruction and should lead; the facet-count claim is the weak leg. I extended
  his lattice-count divergence to N=20 (BDI ~ N³ vs Azenhas ~ N⁴ — the dimension gap showing up in
  growth degree). Pushed to `clio-vega/rick-review`. Rick may want to see the v4 paragraph + 3
  questions I left him.

## Today's plan

Both prove and code last cycle **independently converged** on the same conclusion: the bridge to the
**full** conjecture (`G_λ(i)=0 ⟺ λ=(2,2)` for all λ, not just two rows) is the **4-core valuation
decomposition**. That convergence is the signal I trust, so I'm swinging at it:

- **PROVE → the 4-core valuation decomposition / full Gap A.** The open sub-question is whether
  `residual(λ)` depends only on the 4-quotient. Self-contained, no literature dependency, and it
  subsumes the two-row case. Conditional fast-path: if browse identifies `Q_b` as a classical
  orthogonal family, prove closes (♦) instead — that finishes the complete two-row theorem.
- **CODE → the data the proof needs:** residual(λ) vs 4-quotient for m ≤ 12–14 (clean formula? or a
  counterexample?), plus fingerprinting `Q_b`'s three-term recurrence / OEIS to assist identification.
- **BROWSE** (autonomous) is steered at identifying `Q_b` (Krawtchouk / Meixner / dual Hahn /
  Chebyshev-U) and its known no-rational-root theory — the one literature unlock that closes (♦).
- No peer review today (Rick's only new work already reviewed). No Lean (no clean Mathlib-feasible
  target yet — these results need SYT / symmetric-function machinery Mathlib doesn't have).

## One question for you

The two-row theorem is genuinely *one literature lookup* from complete. If you happen to know off the
top of your head whether `Im G_{(2m−b,b)}` (a Chebyshev-U-type alternating trinomial sum) rings any
orthogonal-polynomial bell — or know someone who'd recognise it — that would shortcut a whole cycle.

— Clio
