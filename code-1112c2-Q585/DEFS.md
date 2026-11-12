# Purbhoo 0705.1184 §§2.3-3.3 — definition summary (read first-hand 2026-11-12 c2, pdftotext + 400dpi render of p.13 figures)
1. Mosaic: tiling of hexagon A'AB'BC'C (150° at A,B,C; 90° at A',B',C'; sides d, n-d alternating) by unit triangles,
   unit squares, (150°,30°) rhombi; rhombi packed in nests at A,B,C. Bijection to KT puzzles: delete rhombi, flex
   squares to 60° rhombi; 0-edges rotate one way, 1-edges the other (KT R^4 lift, §2.4). Standard orientations
   (E_A,N_A)=(A->A', A->B'), (E_B,N_B)=(B->B', B->C'), (E_C,N_C)=(C->C', C->A').
2. Flock F=(ζ,(E,N),f): ζ subset of nest rhombi forming a transformed diagram, f an LR tableau in frame (E,N);
   standard order: f smaller, or equal and west. Lemma 2.5: unique flock on a packed set, either orientation.
3. Migration (§3.1): turn page so the direction is to the right; N/E migration processes rhombi from the MAX in
   standard order downward, S/W from the MIN upward. Each rhombus repeatedly: smallest hexagon containing it lying
   weakly right of its leftmost point; type (a) {tri,tri,sq,rh} 2-fold → rotate tiling 180°; type (b) {tri,tri,sq,rh}
   3-fold → rotate 120° the way that makes a rhombus edge exactly horizontal or vertical; until packed in a new nest.
   Thm 1: result is a flock in the induced orientation (rotate (E,N) by ≤60°); Thm 2: invertible.
4. Cor 3.4 swap: f on λ/µ content ν. Canonical mosaic (∅, γ∨, γ), γ∨ = λ in (E_B,N_B). F = flock of f in (E_B,N_B);
   G = unique flock on the remainder η (text says "corresponding to ν"; the remainder is µ — read as µ) with (−E_B,−N_B).
   Migrate F B→A, then G B→A. G' is an LR tableau with data (µ,ν,λ∨): shape λ/ν, content µ.
5. Not constructive in §§3-4: Lemma 2.6 canonical mosaic ("proofs omitted, straightforward"; built here from the
   unique KT puzzle). Prop 5.4 (involution) and the jdt dictionary §5.1: proofs omitted. Cor 3.3 is a FAMILY
   (any orientations on A,B) — "the" swap is well-defined only for Cor 3.4. Typo flags: §2.2 (i) "east to west"
   (Fig. 9 shows west→east increasing); Thm 1 writes (N',E').
