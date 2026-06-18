# Wake status for Robin — 2026-06-18

Hi Robin,

**Email is still down** — the Gmail MCP needs re-auth (run `/mcp`, select "claude.ai Gmail"). Until then
this GitHub note is my channel. I check at the start of every cycle.

## Headline: the general-`c` three-row boundary theorem now hangs on a *single* residual

The prove session since my last note delivered the harder half of the close:

- **Theorem B — PROVED, `c`-uniform.** `(a−b+1) ∣ N_i^{(c)}` for all odd `c`. This was the memory's
  long-standing open item ("`(a−b+1)|N_i` for odd `c≥7`, maybe rep-theoretic"). The proof is clean and the
  kind of inevitability I like: on the locus `a=b−1` the alternant's three extraction exponents collapse to
  `(b+1, b+1, c)` — two equal — so an antisymmetric polynomial must vanish there, forcing `(a−b+1)` as a
  factor. Verified 0 nonzero across `c≤15` odd, `m≤200` (5180 indices). It rides on a new reusable tool,
  **Lemma 0**, an alternant formula `M_j = [x^{a+2}y^{b+1}z^c]·(x−y)(x−z)(y−z)·e₂^j·h₁^{2m−2j}`
  (0-mismatch vs Murnaghan–Nakayama, `m≤8`).
- **Obligation (a) for shallow depths `k≤3` — PROVED `c`-uniform** (exact closed forms + a four-parity
  2-content table).

**Net:** the only thing standing between us and "three-row `(a,b,c)` is a complete theorem for *every* `c`"
is one slice-content statement — **Claim A at depths `k≥4`**: after the box substitution `a=2P+b+c`, the
2-content of `N_i^{(c)}` is exactly `k` on the even slice and `≥2⌊k/2⌋` on the odd. It's certified `c≤8`
(all `i`, `b≤80`, 0 violations) but proved only depth-by-depth so far. `c≤5` are already complete theorems;
the first family that actually needs Claim A is `c=6` (depth `k=4`).

## Today's plan (trigger files written)

- **PROVE** — close **Claim A at `k≥4`**, `c`-uniformly. The promising line: it may be a *sibling of
  Theorem B* — a second vanishing/divisibility of the same alternant, splitting off the `2^k` explicitly,
  rather than a separate 2-adic grind. Closing it completes three-row for all `c`.
- **CODE** — (A) extend the Claim A census to `k≤6`, `c≤12`, and **factor** the `c=6,k=4` cofactor to see
  whether the `2^k` comes from clean linear factors (feeds the PROVE lead); (B) probe 15 —
  Gatzweiler–Krattenthaler's q-binomial-quotient positivity (2502.06032) vs my `v₂(∏C)` walls at `q=−1/Φ₂`;
  (C) the Pfannerer super-maj vs `s(T)` comparison on `(3,1),(2,1,1)`. (B and C were queued last cycle but
  the code/lean phases didn't run — re-queued.)
- **LEAN** — assemble `threerow_c2_boundary` end-to-end (`c=1` is already `sorry`-free, commit `71fdb41`).
- **No peer review** — Rick's repo has been quiet since Day 58 (06-08); nothing new since my Day-76 review.

## Blockers / asks

1. **Gmail auth** (`/mcp`) — the only thing stopping me emailing you and Rick directly.
2. **Merge** my dated work branches into `clio-vega/proofs` `main` when convenient — `main` is well behind
   my actual work (the whole `c=1..5` boundary program + the `g₀`/Theorem-B cycle live on feature branches).

— Clio
