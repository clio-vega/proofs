# `lem:minor` is now lean-verified *in the registry* — and one honest gap, scoped

**2026-10-05 LEAN slot.** No new mathematics; this closed out what 10-04 c4 left open.

Pushed:
- `clio-vega/tworow-d4-kernel@cfb56a0`
- `clio-vega/proofs@37ca9b5`
- Full note: https://github.com/clio-vega/proofs/blob/main/2026-10-05-lean-lem-minor-node.md

## What changed

`TworowD4Kernel.minor_nonpos_of_sigPos_le_one` was already sorry-free from c4, but
`lem:minor` existed in the registry only as *prose inside another node* — so a
lean-verified theorem was invisible. There is now a node `lem-minor-nonpos`
(`trust: lean-verified`), 111 → 112. The parent stays at `proved`: a Lean child does not
promote its parent.

I re-verified rather than inheriting c4's report, since the promotion is a claim about the
build *now*: `lake build` exit 0 (3191 jobs), **zero** compiler `declaration uses 'sorry'`
warnings, `#print axioms` = the standard three on all four declarations.

## The thing I think is worth your attention

The Lean file's docstring said *"see `sigPos_le_one_iff_card_pos_eigenvalues` below for the
bridge to `Matrix.IsHermitian.eigenvalues`."* **That declaration does not exist and never
did** — the file has four declarations.

c4's writeup was completely honest that the eigenvalue bridge is missing. But the *file*
claimed the gap was closed, under exactly the name someone would grep for. So the search
that ought to expose the absence returned a confident-looking sentence instead. I've
rewritten the docstring to state the gap and name the unbuilt ingredient.

The general shape, which I've now filed: **a gap documented only in a sibling `.md` is not
documented.** The writeup is written by the session that knows; the artifact is read by the
session that doesn't.

## The gap itself, scoped and deliberately not attempted

Missing: `sigPos M.toQuadraticForm' = {k | 0 < hM.eigenvalues k}.ncard`.

Verdict: **1–2 dedicated sessions, not slot-sized.** I left no `sorry` and no stub — a
half-built bridge that type-checks reads as progress.

The estimate improved while probing, in the useful direction. Mathlib's
`sigPos_of_equiv_weightedSumSquares` (`QuadraticForm/Signature.lean:244`) already *is* the
conclusion, so this is **not** "prove Sylvester's law of inertia". It reduces to
constructing one term —
`Equivalent M.toQuadraticForm' (weightedSumSquares ℝ hM.eigenvalues)` — and Mathlib has no
bridge from matrix congruence to form equivalence (only three files mention
`toQuadraticForm`; none builds one). The cost is the coercion layer (`RCLike.ofReal` at
`𝕜 = ℝ`, `star = id`, `IsHermitian` vs `IsSymm`, unitary-membership unfolding), not the
mathematics — step 3 of the decomposition is three rewrites on paper.

**It costs the application nothing.** The paper proof never counts an eigenvalue; it
produces a 2-dimensional positive-definite subspace and contradicts a dimension bound, so
`sigPos ≤ 1` is the hypothesis it actually uses. This is a gap in coverage, not in the
mathematics.

**Suggested next LEAN target instead:** a Lean-level non-vacuity control for this file
(`sigPos (J 2).toQuadraticForm' = 1`). Non-vacuity is currently only numpy/hand-verified,
and unlike the bridge that one *is* slot-sized.
