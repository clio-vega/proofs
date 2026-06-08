# Two-row d=4 order law — code session, 2026-06-08

Computational scaffolding for the two open pieces of the two-row `d=4` order law
`(★)`: the non-vanishing of `I_b(m) = Im[u^b]((1−u)(1+(1+i)u+u²)^m)` for `m ≥ b`.
Reduction: `I_b = (m)_{R+1}·Q_b`, `R = ⌊(b−1)/2⌋`; `q_b = primitive(big factor of Q_b)`;
the law `⟺ q_b` has no rational root.

## Shared module
- `dfour_tworow.py` — builds `I_b, Q_b, q_b`; **independently cross-checked** against a
  direct `u`-expansion. Run `python3 dfour_tworow.py` for the self-test.

## Job 1 — `b ≡ 0 (mod 4)` (supports the PRIMARY prove target)
- `job1_mod4_reduction.py` — the mod-4 reduction is **exact** for `b ∈ {8,…,48}`, and the
  mechanism (`Im((1+i)^j) ≡ 0 mod 4` for all `j ∉ {1,2,3}`) is **`b`-uniform**.
- `job1_fp2_recurrence.py` — the headline 𝔽₂ result: `s=1+i` is **nilpotent mod 2**
  (`s²=2i≡0`), so the `m`-th power truncates, giving the closed form
  `I_b(m) ≡ m·C(m−1,R) (mod 2)` uniformly in `b`. Confirms `q_b ≡ (m²+m+1)^{⌊(b−1)/4⌋}`
  (`b≡0,1`) resp. `m·(m²+m+1)^{…}` (`b≡2,3`).
- `job1_v2_digits.py` — the valuation proposition `v₂(I_b)=v₂((m)_{R+1})−v₂(R!)` holds with
  **no exceptions** over `m ∈ [b,b+400]`; the valuation-minimal `j`-set is **uniformly**
  `{1}` (`m` even) / `{1,2,3}` (`m` odd) — the "3 odds = odd" 2-adic-unit proof skeleton.
- Findings: `FINDINGS-fp2-recurrence.md`, `FINDINGS-job1-b0mod4.md`.

## Job 2 — `b ≡ 2,3 (mod 4)` (supports the SECONDARY prove target)
- `job2_oddprime_newton.py` — coeff factorisations, per-`b` "no root mod p" primes, `disc`.
  (Bounded trial-division, never hangs. Its `all-|slope|<1` flag is **superseded** — see below.)
- `job2_newton_recheck.py` — the **correct** Newton test (all slopes **non-integral**):
  **zero** certificates for monic `q_b`, reciprocal, or scaled — the Newton route is dead.
- `job2_covering_periodicity.py` — covering-set / `b`-periodicity hunt: **no** uniform prime,
  **no** finite covering set (apparent covers evaporate by `b ≤ 62`), **no** `b`-periodicity.
- Findings: `FINDINGS-oddprime-newton.md`.

## One-line verdicts
- `b ≡ 0 (mod 4)`: essentially closed, modulo a finite uniform binomial-valuation lemma.
- `b ≡ 2,3 (mod 4)`: genuinely **non-local**; needs a global (analytic / irreducibility)
  argument. Each individual `b` has a "no root mod p" witness, but nothing uniform exists.

All arithmetic is exact (ℤ / 𝔽₂ / 𝔽_p). No floats, no Sage.
