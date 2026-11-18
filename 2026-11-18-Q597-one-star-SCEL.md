# PROVE 2026-11-18: Q597, one-star SCEL (journey locality)

Full write-up: `2026-11-18-Q597-one-star-SCEL.tex` (+pdf). Code: `code-1118-Q597/`.

**Theorem targeted.** One-star SCEL: (L) Δ_{Mx}(y) = Δ_M(y) off star(u); (E) Δ_{Mx}(y) = s − Δ_M(x) on star(u).
**Status: NOT proved.** Reduced to a journey-level block-swap hypothesis (BS). The reduction is proved.

## Results
- **Step 0.** Cap 60 was never binding: caps 60 and 300 both visit all 5194 states, so the 11-17 data is the *whole* reachable set of the |λ|,|μ|≤5 gl4 census. SCEL with frozen signatures: 1558/1558. A size-6 census run was launched (see the tex addendum).
- **Step 1.** `Board.step` = Purbhoo §3.1 step. Window: tiles within 3 adjacencies of the rhombus; writes 4 tiles. Each step changes exactly 2 pure cells, keeping keys (31,809 steps). Journey delta = Σ step deltas (4674/4674).
- **Computed: re-synchronisation.** On every B/C pair both journeys differ from solo only on an equal-length middle (2–9 steps), and then the rhombus is in the same tile. B/B non-commuting: never (0/440 per journey).
- **Computed: BS.** On the clean journey (y for q-side, x for p-side) the middle is a·C in M and (C−α)·(a+γ) in the other state, with every step a rigid Z⁴-translate. 521/521; B/B 0/880. α = (0,−1;0,0) on q, (0,0;0,1) on p. This is the rhombus 2-cell at the level of journeys.
- **Proved (Thm red, Cor side).** Lemmas window/covariance/additivity plus the prefix–suffix proposition. Under BS: Y^x − Y = Σ_C(T_{−α}δ − δ) + T_γδ_a − δ_a. It is local, and it is one-sided when α is a single-block vector, which explains "one side only". Formula = observed 521/521.
- **Refuted route.** Journeys do not translate their wakes (1520/3116 fail).

## Gaps
G1 BS itself (a crossing analysis à la Purbhoo's Wake Crossing Lemma, for journeys to different nests); G2 (E) needs window data (the middles do not read all of I); G3 the other journey has no translate structure, but thm:cond needs both; G4 catalogue completeness (2 q + 4 p one-star configurations); G5 intermediate boards in 𝓕.
