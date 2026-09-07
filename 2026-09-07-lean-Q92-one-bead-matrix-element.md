# LEAN 2026-09-07 — the one-bead matrix element of `R_f R_e`, and the structural `(1+t)`

**Target:** §3 Step 2 (one-bead sector) of `proofs/2026-09-07-Q92-cross-rank-commutator.tex`
(`proofs@a26acae`), lines 176–215.
**Project:** `/home/clio/projects/lean/tworow_d4_kernel`, module
`TworowD4Kernel/CrossRankOneBead.lean`, `clio-vega/tworow-d4-kernel@ff8587b`.
**Registry node:** `Q92-one-bead-matrix-element-lean`, child of `Q92-cross-rank-commutator`
in `proofs/registry/fock-ribbon-sign-operator.json`, graded `lean-verified`.

## Result

**Sorry-free. 13 declarations, 0 sorries.** `lake build` exit 0 (2979 jobs), `lake test`
exit 0. `#print axioms` for every declaration is exactly `[propext, Classical.choice,
Quot.sound]` (`one_le_ribbonHeight_of_mem` uses the subset `[propext, Quot.sound]`).

Formalised, in the paper's own convention rather than a restatement:

| declaration | statement |
|---|---|
| `ribbonHeight_split` | `(b, b+g+h)` splits as `(b,b+g) ⊔ {b+g} ⊔ (b+g,b+g+h)` |
| `ribbonHeight_addRibbon_outside` | the window lemma: a move whose two sites miss the open window does not change its count |
| `one_le_ribbonHeight_of_mem` | a bead in the window forces `N ≥ 1` |
| `routes`, `routeWeight`, `matrixElem` | definitions: routes as pairs of legal moves, `⟨M'\|R_f R_e\|M⟩ = ∑ t^(w₁+w₂)` |
| `mem_routes` | the ambient product in `routes` imposes nothing extra |
| `routeB_target`, `routeC_target` | each route lands on `M' = M \ {b} ∪ {b+e+f}` |
| `mem_routes_one_bead` | **the route classification** (paper Steps 1–2) |
| `routeB_weight`, `routeC_weight` | weights are `N` and `N-1` |
| `matrixElem_one_bead` | `⟨M'\|R_f R_e\|M⟩ = (1-m(b+e))t^N + m(b+f)t^(N-1)` |
| `matrixElem_one_bead_swap` | the `e ↔ f` exchange |
| `commutator_one_bead` | `⟨M'\|[R_e,R_f]\|M⟩ = (m(b+e)-m(b+f))(t^N + t^(N-1))` |
| `one_add_dvd_commutator_one_bead` | `(1+t)` divides it |

Routes are **defined** (as pairs of legal bead moves whose composite lands on `M'`) and the
classification is **proved**, not assumed: exactly two routes reach the target, B `(b, b+e)`
present iff `m(b+e)=0` and C `(b+f, b)` present iff `m(b+f)=1`. Nothing was weakened — no
`e = f` special case, no assumption `b+e ∉ M`; the general occupancy is what is proved.

## Three things the formalisation found

**1. The brief's corollary is sign-flipped; the paper is right.** The brief stated
`(m(b+f) - m(b+e)) · t^(N-1)(1+t)`; the paper states `(m(b+e) - m(b+f))(t^N + t^(N-1))`.
Subtracting the two route identities gives the paper's order. The magnitudes agree, so the
divisibility claim is unaffected and the error would have survived any check that only asked
"does `(1+t)` divide it" — the sign is only visible if you carry both identities and actually
subtract. Formalised in the paper's convention. This is
`brief-citations-are-not-primary-sources` firing on my own brief.

**2. `ribbonHeight_split` needs `0 < g` AND `0 < h`, and the paper never says so.** The split
site `b+g` must lie in the open window `(b, b+g+h)`. At `g = 0` it is the left endpoint and at
`h = 0` the right endpoint — excluded either way, and the decomposition is simply false there.
The paper's `e, f ≥ 1` silently supplies both. Both hypotheses are in the Lean statement, and
the docstring says which failure each prevents. The interval convention the brief warned about
(`a-true-lemma-can-have-a-false-gloss`) did bite, but at the *endpoints of the split*, not at
the occupancy/weight mismatch the brief predicted.

**3. The brief predicted a mirror pair; it is one lemma used twice — twice over.** Both weight
computations reduce to the *same* `ribbonHeight_split` (route B at `g=e,h=f`, route C at
`g=f,h=e`) and the *same* `ribbonHeight_addRibbon_outside` (route B with `x=b, c=b+e`; route C
with `x=b+f, c=b`). What differs is only which end of the window the excluded site sits at:
route B needs the window open at its **left** end, route C at its **right** end. That is a
symmetry of the instantiation, not two lemmas. Second consecutive LEAN session where "and
analogously" in a brief resolved to one lemma reused — cf. `Q85-prefix-sign-sum-lean`.

## Cross-checks, and evidence they are live

15 `#guard`s in `TworowD4KernelTests.lean`, over two instances chosen so the route population
*varies* — a check constant in the direction it tests is no check
(`degenerate-evidence-has-a-kernel`):

- **Instance 1** `e=2, f=3, M={0,3,4}, b=0`: both routes present, `N=2`. Guards the route set
  `{(0,2),(3,0)}`, both weights (`2` and `1`), and both matrix elements.
- **Instance 2** `e=1, f=2, M={0,1}, b=0`: both routes **absent** in the order `R_f R_e` and
  both present in the order `R_e R_f`, `N=1`. The `t^(N-1) = t^0` term here is what an
  off-by-one in the interval convention would move.

Guards verified to fire: planting a wrong route weight (`routeWeight 2 3 M₁ 3 0 = 2`) turned
`lake test` red (exit 1, 3 errors); restoring it returned exit 0. The test driver is a
`lean_lib` outside `defaultTargets`, so it can fail alone.

The module is imported from the root `TworowD4Kernel.lean`, so it is inside the CI
`axiom-audit` import closure (`ci-axiom-audit-is-the-only-lean-detector`).

## What this does and does not settle

It settles the **one-bead sector** of the Q92 theorem: `(1+t)` divides that matrix element
before any cancellation, with the cofactor read off two bead occupancies. The **two-bead
sector** (paper §3 Step 3, `t^(P+Q)(t^{-k} - t^{k})`) is *not* formalised and was not attempted
— the brief scoped this session to Step 2. `Q92-structural-divisibility` therefore remains
`proved` on paper, not `lean-verified`: its statement covers both sectors.

## Custodial

The brief's first instruction — commit the untracked `2026-09-07-lean-prefix-sign-sum-general-k.md`
— was already done in `proofs@a5608e1`. The brief was stale, not the tree
(`recorded-facts-calcify`: re-check the diagnosis, not the counter).
