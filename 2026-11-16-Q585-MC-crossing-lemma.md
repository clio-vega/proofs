# Q585 (MC) crossing lemma, PROVE 2026-11-16: NOT proved (companion to the .tex)

Full write-up: `2026-11-16-Q585-MC-crossing-lemma.tex` / `.pdf` (commit 81aedd9). Code: `code-1116-Q585-MC/` @ 98f5c79.

- **Theorem targeted:** (MC) x→B;y→C = y→C;x→B (same tiling, same landings). Also (SW). Result: **neither proved**. G1 stays conditional.
- **Proved:** Lemma (degree) and Lemma (hexagon classification). Every 180° step slides a unit lozenge of two same-type triangles by a unit
  vector of the other lattice, and every 120° step is Δ₀+Δ₁. Corollary: the triangle shapes are conserved, so the state's triangles are
  (q,p). Caveat: two planar-topology steps are argued informally. They are corroborated by 6982 + 2420 steps.
- **Computed:** the no-repeat and direction laws (each migration is a unit push of a barycentric height). The union law holds on B/C
  1558/1558 and on SW 987/987, and fails on all 440 non-commuting B/B pairs (control). RIGID' holds 9982/9982. RIGID fails 2/5192.
- **Gap:** a local exchange rule for the triangles in W_B ∩ W_C that keeps the per-triangle sums equal along the crossing window, plus a proof of RIGID' or of the landing half.
- **Dead ends:** transporting the Wake Crossing Lemma (it is a same-direction lemma); the shelf (1209.1262, 1401.1725).
