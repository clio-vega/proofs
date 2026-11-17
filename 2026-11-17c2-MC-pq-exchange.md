# 2026-11-17 c2 — (MC) on pure-cell data: the star-chain exchange law

Full artefact: `2026-11-17c2-MC-pq-exchange.tex` (pdf alongside). Code: `code-1117c2-MCpq/`.

## Statements
- **SCEL (computed):** for an overlapping B/C pair (x→B, y→C) at M, let X, Y be the solo net (p,q)-deltas and Y^x, X^y the deltas of the second move.
  Then there is a segment U of interior lattice points on ONE side, running in direction (1,−1) on q or (1,0) on p, with I = ∪_{u∈U} star(u), such that
  Y^x = Y and X^y = X off I, while Y^x = s − X and X^y = s − Y on I. Here s = b + c is the step shared by B and C. U is recognised from the solo signatures alone, by 7 signatures (single/first/inner/last).
- **Theorem (proved):** SCEL ⇒ (MC), by additivity of net deltas plus PQ-RIGID.
- **Picture:** the hexagon v+star(u) goes from v to v+s across the diagonal of the rhombus {v, v+b, v+s, v+c}. The two orders go round its two sides.

## Verification
- 1558 B/C + 779 B/B pairs (census sampler, cap 60). Tiling equality ⟺ (p,q) equality in 2337/2337 pairs.
- Predicted I = observed I in 1558/1558 B/C pairs (1037 empty, 521 nonempty), and in 779/779 held out. Classes |U| = 1,2,3,4 occur 362/128/30/1, and |add| = |rem| = 2|U|.
- Control: the rule predicts no interaction on all 779 B/B pairs, so it is wrong on all 440 non-commuting ones.
- Brief sanity check: class (0,0) agrees trivially in 258/258 pairs as vectors, but in 0/258 as (old, new) assignments.

## Gaps
- SCEL is not proved. The locality of the second journey needs Purbhoo's Board.step (Q597); the static lemmas (lem:ribbon, Lemma W) do not reach it.
- Not explained: why U is a straight segment, and why only one side interacts. Step 4 (reachable set R) was not run.
