# Q594 / Q585 (MC), PROVE 2026-11-16 c2: MC is NOT transported switching (companion to the .tex)

Full write-up: `2026-11-16c2-Q594-MC-switching.tex` / `.pdf`. Code: `code-1116c2-Q594/`.

- **Step 1 (cactus test).** I used the solo triangle wake as the window. Result: a B/B pair commutes ⟺ its two wakes are disjoint (339/339). All 440 non-commuting pairs overlap, and none are nested. The nested cell is empty, so the test was blind to the conjugation relation. The 339 "commuting interacting B/B" pairs therefore do not interact. MC's real content is the **779 overlapping B/C pairs**.
- **Step 2 (sources, first-hand TeX): null.**
  - Buch thm:mutate says two gash paths born at one temporary piece are disjoint. It covers only the trivial cell.
  - BKPT lem:chainorder is a same-direction Pieri swap, and it is vacuous at one step (cor:label012).
  - ACM App. A is a whole-flock dictionary given by example, with moves omitted.
  - So there is no stepwise dictionary, and step 3 was not run.
- **Step 4a: no.** The barycentric height is not an aura. An aura is conserved along a path; the height increases.
- **Step 4b: P-RIGID (computed).** p alone determines the mosaic, and so does q alone. 0 collisions on 1979 mosaics pooled by hexagon and on 3068 pooled by ν. The nest-data control collides 141 and 167 times.
  - **Proved:** φ is injective on Z⁴; the ribbon lemma (no closed ribbons, by a translate-area argument; R_f is a path from p(σ) to p(σ')); no mixed hexagon flip.
  - **Gap:** degenerate fibres p⁻¹(v)=∅.
  - **Conditional:** if P-RIGID holds, then MC ⟺ its 1-triangle half.
- **Exchange data (necessary conditions only).** On overlapping B/C pairs, the second journey drops only first-wake triangles and adds exactly as many as it drops. The "same change" law is false. B/B pairs are unbalanced.
- **Status:** MC, SW and P-RIGID are not proved. G1 stays conditional on (MC)+(SW).
