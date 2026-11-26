# PREDICTION — PROVE 2026-11-26 c2, written before step0s2.py / zone.py were run

Notation (G3prime tex §1): H_E, H_L = the four tiles removed by E_0 = J_E[p_E], L_0 = J_L[p_L] (M0 tiles).
C_E0, C_L0 = tiles created by E_0, L_0. For a step s, G(s) = removed tiles other than its rhombus.
Oth(E_1) := G(E_1) \ (H_E ∪ H_L), Oth(L_1) := G(L_1) \ (H_E ∪ H_L).

## (Share_2), the predicate (per diverging B/C record)
 (S2.i)   G(E_1) ⊂ M0 (no tile created by E_0 other than E_1's rhombus), |G(E_1) ∩ H_L| = 1, G(E_1) ∩ H_E = ∅.
          Symmetrically G(L_1) ⊂ M0, |G(L_1) ∩ H_E| = 1, G(L_1) ∩ H_L = ∅.
 (S2.ii)  Oth(E_1) = Oth(L_1), a set of two M0 tiles.
 (S2.iii) |C(E_1) ∩ C(L_1)| = 2 in case b, 0 in case a.
Prediction: all three hold 779/779 at size 5 (full) and on every record of fresh size-6 chunks 11-12 (of 24).
Prediction on kinds: Oth = {square, triangle}?  I guess one square + one triangle (each G has 1 square, 2 triangles, and
the hex_L tile is a triangle).  Low confidence.
Translation classes of the joint tuple (H_E, H_L, G(E_1), G(L_1)) relative to the shared square: 4 in case b (= W_F
classes), and in case a the number of distinct (E_0,E_1,L_0,L_1) prefixes of the W_a classes (6 at size 5).

## (b) DREAM's zone falsifier
alpha (BS(L) vector) is a lattice unit vector in one factor. Claim Z: every tile removed by a corridor step of E^L or of
L^E has an edge parallel to alpha (two of its points differ by ±alpha), and these tiles form one edge-connected strip
(consecutive ones share an alpha-edge).  Prediction: Z holds on all records with k>0 (256 at size 5).  Confidence ~60%.

## (c) B/B control
On non-commuting B/B pairs (L := role at a class-1 a-step of its target at p, as in g3red.py) S2.ii should FAIL on most
records (prediction: fails on > 50% of the 440), or S2.i is undefined/false.

## Added 18:45, before corr.py was run: the corridor (Merge on M0)
In case b, at every shared rhombus R_i (i>=2) E_i and L_i remove the same tiles. While the type is a, the steps are identical.
The first b-step is (0,b) for both: same removed tiles, different created tiles. Corridor length k == m_L - 3 on every record.
Mechanism: finality puts all removed non-rh tiles in M0; class-0 G_t identical for B and C; Lemma ovl => types agree.
(H2), both journeys step at every shared rhombus: predicted to hold on every record. Case a: E_2 and L_2 at different rhombi.
