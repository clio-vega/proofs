# Q585 re-aimed — Purbhoo's Cor 3.4 flock swap M vs the crystal commutor C (PROVE 2026-11-12 c2)

**Statement (computed, not proved).** For every gl_4 cell with c^ν_{λμ} ≥ 2 and |λ|,|μ| ≤ 5, and every cell with
c ≥ 3 and |λ|,|μ| ≤ 7 (190 cells, 518 elements, c = 2,3,4,5), and for each of two bijections β (β_Q: recording tableau;
β_CM: content-matrix transpose), M ∘ β = β ∘ C on hw(B(λ)⊗B(μ))_ν. M∘M = id (Prop 5.4) with 0 failures.

**Instrument.** An exact migration engine built from Purbhoo's definitions: KT puzzle → mosaic via the R^4 lift,
with Z^4 coordinates. The geometry convention is forced (1 of 4 valid). It was validated against the paper's own
theorems before Cor 3.4 was touched: 139 mosaics × 8 Table-1 variants, Thm 1 holds 100%, and Cor 3.2 is injective.
Flagged interpretation: a tie-break for the 956 non-unique minimal hexagons, read off Figs (a)/(b). It is load-bearing
(the opposite choice is stuck on 63/142).

**Controls.** Lex ≠ C on 190/190 cells. Reverse-lex is near-blind (it differs from C on only 27/123 c≥3 cells). The
c=2 cells carry 1 bit each, because a c=2 cell has only two bijections. Lesion (G with the standard orientation)
fires 67/67.

**Gap.** There is no proof. A proof needs (1) the B→A migrations of Cor 3.4 written as a jdt composite (Purbhoo §5.1
omits the proofs and covers only B→C), and (2) that composite identified with C. The tableau-switching route is recalled
from memory only and is not cited.

Artifact: `proofs/2026-11-12-Q585-purbhoo-swap-vs-commutor.tex` (pdf). Code: `proofs/code-1112c2-Q585/` @ 7303d1db5f52.
