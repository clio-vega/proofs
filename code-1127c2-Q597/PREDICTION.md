# PREDICTION — PROVE 2026-11-27 c2, written BEFORE any run of dict_falsifier.py

Reading done first: Purbhoo 0705.1184 seed pdftotext §2.3 (flocks, accessible), §2.4 (mosaic def), §3.1 (migration, Thm 1, Thm 2);
engine code-1112c2-Q585/migrate.py (Board.step, Board.journey) and census move() (pop nest-A corner, journey on a copy).

Observation from reading (before computing): the engine's `step` IS Purbhoo's per-rhombus operation (minimal rh+sq+2tri hexagon
weakly right of the leftmost point, rotate 180 or 120 with an edge horizontal/vertical), and `journey` stops exactly when the rhombus
is packed into a nest other than its source with the nest still an ideal. A nest-A corner is a removable corner, so {x} = alpha/alpha1
with alpha1 = alpha minus x straight: the singleton flock is ACCESSIBLE by definition. Direction u or w = a side of the nest = E or N
in the frame orientation, which is the "opens NE => north or east" consistency rule. So the dictionary is "J_x = Purbhoo migration of
the accessible singleton flock {x} to nest T" by construction, EXCEPT at ties (registry: 956 ties at n=5 contradict Purbhoo's
uniqueness parenthesis; resolved by our flagged rule).

Also from reading: Thm 1 is NOT needed for "M0.J_x is a Purbhoo mosaic". That is definitional: hexagon rotations re-tile the same region
by the same three shapes, only x moves, and x stops packed. The real content of (Mos) is the converse half of Purbhoo §2.4,
(Surj): every Purbhoo mosaic is M0(ell) for a KT puzzle ell. Thm 1 is about labels (flocks), which a singleton does not carry.

Predictions for the census run (1979-tiling size-5 closure, every journey from every state, both targets):
P1. Start rhombus accessible (removable corner): 100% (by construction of corners()). Confidence 99%.
P2. Target reached = intended nest (B for w, C for u): 100%. 95%.
P3. Every step's rhombus centroid moves weakly right in the migration frame ("roughly to the right, sometimes strictly vertically"): 100%. 85%.
P4. Ties (more than one candidate hexagon) occur on a positive fraction of steps: YES, a few percent. 80%.
    These are the only steps where "Purbhoo migration" is underdetermined; the dictionary holds there only under our tie rule.
P5. Every reached board is a Purbhoo mosaic (tiles cover H exactly: area count + no interior overlap; every nest an ideal): 100%. 97%.
P6. Lesion (pop a NON-corner nest-A cell, run its journey to B): the result is never a mosaic (nest A not an ideal), and bridge-type
    check (labels read off by psi, rebuilt by build_mosaic) fails on >= 90% of lesioned boards. Control (corner pops) passes 100%. 75%.
Violations of P1-P3, P5 would kill the dictionary as stated.
