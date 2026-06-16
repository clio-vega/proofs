# FINDINGS — Job B (2026-06-16 code session): Pfannerer super-maj vs `s(T)`

**One-line verdict.** Both comparisons **FAIL**. Pfannerer super-maj (plain-index sum over
super-descents) **cannot** reproduce the parity-twisted multiset `{s(T)}` — not for any
sign-set rule, and not even for any *per-tableau* choice of `D` (Hall-infeasible). And the
`t=−1` content product `∏_□(1−q^{c(□)})` is **identically zero** for every shape (the diagonal
content-`0` cell contributes `1−q⁰=0`), so it carries no information about `J*`. **Super-maj is
pruned** as a model for `{d_j}={s(T)}`. The leg stays blocked on this route; the next step is the
4th literature home, **not** super-maj.

`s(T) = Σ_{i∈Des(T)} w_i`, `w_i = 2i−1` if `n−i` odd else `0`; `G_λ(q)=Σ_T q^{s(T)}`.
Pfannerer (per CODE.md / 2603.16598): signed tableau `(T⁺∈SYT(λ), D⊆{1..n})`,
super-descent `i` iff `(i∈Des(T⁺) ∧ i+1∉D) ∨ (i∉Des(T⁺) ∧ i∈D)`, `supermaj = Σ` super-descent `i`.

---

## Comparison 1 — super-maj ≠ `s(T)` multiset (decisive NO, all shapes)

Tested `λ ∈ {(3,1),(2,1,1),(2,2),(3,2),(2,2,1),(4,2),(3,1,1),(3,2,1)}`. For **every** shape:

- none of the natural rules `D∈{Des, ∅, all, parity-flip, non-descent∧parity}` matches `{s(T)}`;
- the **Hall feasibility test fails**: even allowing each `T⁺` its own optimal `D`, the achievable
  super-maj values cannot be matched (as a multiset) to `{s(T)}`.

**Why — the mechanism (`λ=(3,1)`, the make-or-break shape).** `{s(T)} = {0,1,5}`. Per tableau:

| `Des(T⁺)` | `s(T)` | achievable super-maj over all `D` |
|---|---|---|
| `{3}` | 5 | `{0,1,2,3,4,5,6}` |
| `{2}` | **0** | `{2,3,4}` |
| `{1}` | 1 | `{1,2,4,5}` |

The value `0` is needed and is reachable **only** by the `Des={3}` tableau. But then the `Des={2}`
tableau must supply a value in `{1,5}` — and its super-maj is confined to `{2,3,4}`. Infeasible.

The structural reason: the **parity twist zeroes even-`(n−i)` descents** (`w_2=0` because `n−2=2`
is even, giving `s=0` for the `Des={2}` tableau), but Pfannerer super-maj **charges that descent
the full index `i=2`** regardless of `D` (its band is `{2,3,4}`, bounded below by 2). No sign-set
can convert a charged index back to `0`. `s(T)` is **sparse and gappy** (`{0,1,5}` — the parity
twist creates the jump `1→5`); super-maj is **dense and unimodal** (`{0:2,1:6,2:10,3:12,4:10,
5:6,6:2}`). Different shapes; no `D` reconciles them.

`(2,1,1)` is the same story: `{s(T)}={1,5,6}`, the `Des={1,2}` tableau has `s=1` but super-maj
band `{2,3,4}` (can't reach 1 alone while the others cover 5,6).

## Comparison 2 — `t=−1` content cancellation is degenerate (NO `J*` model)

`f^λ(q,t)=f^λ(q)∏_□(1+t·q^{c(□)})`; at `t=−1`, `∏_□(1−q^{c(□)})`. **Every partition has a
diagonal cell `(0,0)` of content `0`**, contributing the factor `1−q⁰=0`, so the product is
**identically `0` for every shape** (verified `(3,1),(2,1,1),(2,2),(3,2),(3,2,1)`). The
vanishing order = #content-`0` cells = **Durfee side `d(λ)`**, which does **not** match
`ord_{q=−1}G_λ`:

| `λ` | Durfee side | `ord_{q=−1} G_λ` |
|---|---|---|
| `(3,1)` | 1 | 0 |
| `(2,2)` | 2 | 0 |
| `(3,2)` | 2 | 0 |
| `(3,2,1)` | 2 | **3** |

No alignment. The `t=−1` specialization of the hook-content engine collapses (diagonal kills it)
and gives no tableau model for `J*`.

---

## Verdict & next step

- **Comparison 1: NO** (Hall-infeasible, all shapes). The plain-index super-maj cannot carry the
  parity twist `w_i∈{0,2i−1}`. **Super-maj is finally pruned** for `{d_j}={s(T)}`.
- **Comparison 2: NO.** `t=−1` content product is identically zero; vanishing order = Durfee side
  ≠ `ord_{q=−1}G_λ`. No `J*` tableau model here.
- **The `{d_j}={s(T)}` leg stays blocked on the super-maj route.** Per the dream's own framing,
  the standing block-lifter has now been run and refuted. The next move is the **4th literature
  home** (the triple-convergence list — Pfannerer/Armon–Swanson super-maj is now *removed* from it;
  remaining live candidates: Keith–Anible fixed-descent-count and Prasad eigenvalue-multiset), or a
  direct multiset-equality proof of `{d_j}={s(T)}` that does **not** route through super-maj.
- A genuine model must respect the **gappy, parity-twisted** spectrum of `s(T)` (the `w_i=0` on
  even `n−i`). Any candidate statistic whose values are a *dense* index-sum is excluded by the same
  Hall obstruction.

## Files
`code/pfannerer_supermaj.py` (both comparisons + Hall feasibility + achievable-set diagnostics).
Reuses the `s(T)` machinery from `code/2026-06-13-jobA-oddcontent.py`.
