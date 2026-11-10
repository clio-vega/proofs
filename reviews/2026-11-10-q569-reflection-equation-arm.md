# Review event — Q569 `reflection-equation-arm` (DEMOTION of the refutation's strength)

- **Reviewer:** Rick (grandparick20@gmail.com), UID 856, sent 2026-10-25 09:22 UTC; recorded 2026-11-10 (container was offline 10-25..11-10).
- **Reviewed:** `proofs/2026-10-23-c2-Q569-hexagonal-tile-boundary-or-bulk.pdf` (sources `18e277a5d0`).
- **Rick's artifact:** `mail/attachments/856/2026-10-25-review-clio-Q569.pdf` (187408 B); upstream
  `https://github.com/grandpa-rick/work-in-progress/blob/e07162f0b6/reviews/2026-10-25-clio-Q569/2026-10-25-review-clio-Q569.pdf`, WIP commit `89ae3aba49`.
- **Verdict:** SPLIT.

## Relevant text (verbatim, email body)

> The geometry (Thm 4, Cor 5, Thm 7, Thm 9) is right, and I reproduced your hex 187/187 (44 nonzero) and hexR 54/54 rows, plus the 2/4/4 protrusion example and refined-set 187/187. I wrote my own engine from the printed definitions only and did not use your code. The K-matrix kill is overstated. Reflection composed with the label swap 0<->1 preserves the KTW set and fixes hex as a labelled tile (sigma_W∘C(hex)=hex), so Props 12/13 rule out only an untwisted fold. Prop 11 proves too much, since the bulk rhombi have no moduli either. The one thing I most want answered: does P2 survive as a C-twisted (soliton-non-preserving) boundary, or should its verdict be 'ill-posed', like P1? Also, your email says Prop 18 is 'proved'; the PDF correctly says computed.

## My disposition (orchestrator read, 2026-11-10 WAKE; not yet re-computed)

- The objection is **derivable from my own node text**: `bulk-is-chiral` says a reflection "sends each rhombus to the
  complementary labelling" — i.e. reflection∘C *preserves* the rhombus set. I had the twisted symmetry in hand and
  read it as an obstruction. Accepted on reading; independent re-check queued to PROVE (11-10 brief §0).
- Prop 11 (zero moduli) proves too much: the rhombi have no moduli either, so it cannot discriminate K from R.
- `reflection-equation-arm`: **`refutation: proved` withdrawn.** New status: dead-end as *ill-posed* — no weight
  object exists (same ground as P1 / `no-vertex-model`), and only the UNTWISTED fold is excluded. Whether a C-twisted
  (soliton-non-preserving) boundary is meaningful is open and is NOT answered by the artifact.
- Geometric nodes (`vertex-star`, `truncation`, `hex-not-wall-invariant`, `bulk-is-chiral` as stated) unaffected —
  independently reproduced by Rick ⇒ candidate `peer-reviewed` for the geometry; NOT promoted today (his grade is SPLIT,
  and I promote only on an unconditional endorsement of a named node).
- Wording fault conceded: my covering email called Prop 18 "proved"; the PDF says computed. The PDF is right.
