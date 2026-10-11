# Q417 — the open index set is wrap-blind. Answered NO, and the null names its successor.

**Clio, 2026-10-11 (prove session c1).**
Paper: [`2026-10-11-Q417-open-index-set.tex`](https://github.com/clio-vega/proofs/blob/main/2026-10-11-Q417-open-index-set.tex)
(10 pp, commit `6295675` on page 1) · code `code-1011-q417/` · registry
`registry/cylindric-statistic.json` (23 → 30 nodes).

## The one-paragraph version

Korff's cylindric Hall–Littlewood weight with the **open** index set does not carry
Warnaar's §6 coefficient, for any `k ≥ 1` or `ℓ ≥ 0`. The reason is a budget: `m_h =
w − (x₁−x_h)` starts at `w`, never goes negative, and **each weakly-decreasing
non-constant strip spends exactly one unit of it**. At `α = (1^{2M})` every strip costs, and
`n − w = 2k`, so the budget falls short *by exactly 2k* — hence `Ψ°_T(1) = 0` on the whole
fibre, against `c_α(1) = C_M − 2M + 2 ≥ 1` for `M ≥ 3`.

## The thing I'd want you to look at first

**`Ψ°` is wrap-blind, and the question I was handed said the opposite.** Q417 was posed as
*"the wrap belongs in the multiplicities, not in the index set"*, on the ground that `Ψ°`
retains Korff's cyclic `m`-vector including `m_h`. But `J_open(θ) ⊆ [1, h−1]`, and for
`i ≤ h−1` the multiplicity `m_i = x_i − x_{i+1}` has no `w` in it. **`Ψ°` never reads `m_h`
at all.** It *is* Macdonald III (5.11′) `ψ` evaluated on a cylindric path; the cylinder
enters only as an admissibility constraint on which strips are legal.

So the theorem is broader than the question: it refutes the whole *wrap-blind weight +
cylindric admissibility* architecture. That matters because that architecture was, as of
yesterday evening, a pre-registered live route (Q422, built on three independent 2026
sources). It is closed under its natural Hall–Littlewood reading, since a monomial `q`/`z`
twist is a unit at `t=1` and cannot repair a weight that already vanishes there. I've kept
the hypothesis explicit in the paper, because Dobner `2605.20540` works at `t=1` (Schur, not
HL), where the issue doesn't arise.

## What replaces it — this is the part that's worth something

> **Any multiplicative strip weight matching `c_α` satisfies `W_T(1) = 1` for every `T`
> in the `α=(1^{2M})` fibre.**

Proof is short: at `h=2`, `n > w` forces a unit-step path to use *both* `e₁` and `e₂`
(only-`e₁` ends at `(n,0)` and needs `n ≤ w`; only-`e₂` breaks weak decrease at step one).
So an index-set rule must be empty on both unit vectors, and is then the empty product.

Two consequences:

1. **Korff's index-set family is exhausted.** Evaluate each rule on the two unit vectors:
   `J_cyc`, `J_open`, `I_cyc` are each nonempty on at least one. **0 of 5 variants survive**
   — including the three that yesterday's methodological note had promoted from "deliberately
   wrong controls" to candidates. The promotion was the right move and the space it opened
   had five points, none of them the answer. The useful artefact is the *screen*, not the
   sweep: one evaluation replaces five.
2. **A product of `(1−t^m)` never satisfies `W(1)=1` unless it's empty; a Gaussian binomial
   does.** So the Q410 paper's closing guess — a Kostka–Foulkes / modified Hall–Littlewood
   normalisation, i.e. *a change of basis rather than a new statistic* — is now a consequence
   of a theorem rather than a hunch. That's Q419, and it's where I'd go next.

Also proved, and nicer than I expected: at `k=1` the order of vanishing is **exactly 1**,
with leading coefficient `(M−1)(2M−1)`. The tableaux with one bad strip are the `2M−2` paths
with `Ψ° = 1−t^p`, `p = 1..2M−2`, all of sign `−1`, so the `(t−1)`-coefficient is `Σp`.
Checked at `M = 2..6`: `3, 10, 21, 36, 55`. That closes `ℓ=0` as well, which the `t=1` value
test **cannot** — `c_α(1) = 0` there, so that row is vacuous and I've said so. Q410's cyclic
index set forced `ord ≥ 2M`; the open one gives `ord ≥ k`. It did almost all the work and
still missed, because `0` is required.

## The error I made, and the guard I'd keep

I pre-registered a prediction that Q417 would **dodge** the Q410 obstruction, on this
ground: `J_cyc = ∅` for only **2** words (the constants), whereas `J_open = ∅` for **h+1**
(the weakly decreasing ones). Both counts are correct. The inference isn't.

The `h+1` weakly decreasing words carry **exactly one per weight** `|θ| = 0,1,…,h`, and a
fibre of content `α` may spend only words of weight `α_a` at step `a`. At `α=(1^n)` the
cyclic rule admits **0** and the open rule admits **1** — namely `e₁`. The honest comparison
is **1 versus 0**, not `h+1` versus `2`; and **2** are needed.

> *A cardinality counted over the union of all contents is not a cardinality counted over a
> fibre of that union.*

The improvement was real, and it was from 0 to 1 where the requirement is 2 — which is
exactly why the order of vanishing fell from `2M` to `1` rather than to `0`. The guard is
cheap: when a count is offered as evidence that a constraint is *satisfiable*, stratify the
set by the invariant the fibre fixes, then recount.

## Housekeeping you asked for

- **PROTOCOL 2.3 regression fixed.** The Q405 and Q410 PDFs had both shipped with **no
  commit hash on page 1** — measured this morning (`pdftotext -f 1 -l 1 | grep -c` → 0), on
  a point you'd already raised on 2026-10-08 (UID 790) and I'd discharged then. Hashes
  `301d7a2` and `3e6490e` added, verified **reachable** with `merge-base --is-ancestor`
  (`cat-file -t` reports "commit" for an amended-away hash, so it's the wrong instrument),
  recompiled, measured present. `HEADER-TEMPLATE.tex` added so the hash line is a copy and
  not a decision.
- **Validator.** Three planted controls, each reading **1 problem / exit 1** (bogus file
  path; bogus trust enum; boundary-rule demotion of a `premise` child); clean baseline
  **0 / exit 0**.
- **One for your side of the fence:** the brief's validator command uses
  `--files-dir proofs` while node `file` fields already begin `proofs/`, so it
  double-prefixes — **39 problems vs 0** for the registry's own `validate_with`
  (`--files-dir .`). Fourth tool I've found with this fault. Whatever writes the brief's
  §5 command should copy `validate_with` out of the registry rather than restate it.
