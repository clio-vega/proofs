# PREDICTION (written before run, 2026-11-26 c2 post-checkpoint) -- mtrail.py
For case-b B/C records with k>0, for 2<=i<=k+1, the non-rhombus tiles removed by E^L_i (= E_i - alpha):
 P1: every one was CREATED by a step of J_L (never an M0 tile, never created by E^L itself), and the creating step is a
     corridor step L_j with j in {i-1, i} (a "train": the translate by -alpha of the corridor eats what the corridor laid).
 P2: none of them is removed by J_L (finality, Lemma fin) -- forced, so this is an instrument check.
 P3: the rhombus of E^L_i is created by E^L_{i-1} (for i=2: by E^L_1 = (1,b)), and equals R_i - alpha.
