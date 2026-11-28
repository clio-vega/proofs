# Q597: journeys terminate — (NS!) proved (PROVE 2026-11-28 c1)

Full write-up: `2026-11-28-Q597-termination.tex` (8 pp, compiled). Code: `code-1128-Q597/`.

## Theorem
Let M be a Z^4-labelled Purbhoo mosaic (class M, 0<m0<n). Pop a nest-A corner and run its journey toward T in {B,C}. Then
at every state that is not admissibly packed in a nest other than A, **exactly one candidate is good and every other one
is bad**. So (S2) leaves the good one and (S3) applies it. Hence every journey terminates within
#(non-rhombus tiles)/3 steps (Prop term, 11-27 c2), and **(MC) holds for every B/C pair, unconditionally**.

## Proof (outline; every step is in the tex)
1. Brief's sharpening: goodness is free, because (S2) always drops a bad candidate when there are >=2. Uniqueness is NOT free.
2. Shadow lemma: on an intermediate board, F = psi o phi^{-1} (Thm S) still satisfies (a) nest rhombi -> sides,
   (b) dH -> dT_n bijectively, (c) F(H) in T_n (the degree count never used where rhombi sit), and (d) R -> the gash [P,Q],
   with both obtuse corners -> M.
3. Lemma cand (Euler, V=7, one interior vertex): every candidate is R + the full fan at an obtuse corner.
4. Lemma fan (finite, 72 configurations, the engine's own code): 210 = 60+60+90 gives 3 fillings. Exactly one corner v*
   is good for all 3 fillings; the other is bad or fails the x-rule; v* is the obtuse end of the across edge.
5. Case I, M interior: both fans are triangles and squares, so exactly one good candidate.
   Case II, M on dT_n: the gash lies along side s; the outer edges are the back edges of an addable frame cell, so R is
   packed in s. If s = A, the single candidate is at v_in = v*, so it is good.
6. Induction with Lemma class (W1).

## Verification
Size 5 (5194 states, 65802 pre-step states): exactly 1 good candidate 65802/65802; bad = undo 19207/19207; the engine
removes the good one 65802/65802; Case I 56060, Case II 9742 (all k=0); T_n tiling 65802/65802.
Size 6, fresh chunk 19: 21353/21353 the same.
Lesions: L1 non-corner start stuck 1638/1638 (as predicted); L2 a deleted tile gives 359 stuck + 147 cycling of 1984;
L4 x-rule off gives 2 configurations at Case I states and 1 at Case II, exactly as predicted.

## Gaps
None open. Computer-assisted finite lemmas: Lemma fan (ii)-(iv), (Hex) (11-27 c2), Lemma table(2) (11-20).
