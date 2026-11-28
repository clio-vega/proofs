# PREDICTION — PROVE 2026-11-28 c1, written BEFORE any run (ns_census.py)

Reading done first: no-repeat tex (Lemma class, table, undo, start, N, Prop good), MC tex (Def M, Lemma mosP, Thm S,
Prop term), engine migrate.py (Board.step / journey).

## Reading of the brief's sharpening (before computing)
(NS) <=> "exactly one GOOD candidate at every non-packed state reached by good steps".  NOT pure existence:
two good candidates both survive (S2) (they contain the same across tile), and for the class with undefined across
edge (S2) does nothing; either way len(cands)>1 -> RuntimeError (stuck).  Bad candidates are harmless only because
(S2) removes them when a good one exists; if the across tile does not exist and a bad one is the only... (bad has
across edge on dH, so it never contains the across tile; S2 drops it whenever |cands|>=2).
So the brief's "(NS) is pure existence" is half right: goodness is free, uniqueness is NOT free.

## Idea to test (route 1, but corrected)
Intermediate boards are NOT in M (the moving rhombus R is not in a nest).  But R = s0 x s1 with s0, s1 the same
vector u (Lemma class), so psi collapses R to a length-2 segment [P,Q] = [psi b, psi b + 2u] ("gash"), and the
degree argument of Thm S Steps 1-4 never used where R sits (a rhombus contributes 0).  Prediction (G):
G1. tri/sq images tile T_n exactly (each unit triangle once): 100% of pre-step states.  Confidence 90%.
G2. labels agree across every interior puzzle edge except the 4 half-edges of the gash; along the gash the two
    sides read (0,1) and (1,0).  100%.  Confidence 80%.
Then (NS) becomes a statement about a KT puzzle with a gash, which is LOCAL in the puzzle picture.

## Census predictions (size-5 uncapped closure, every state, every nest-A corner, B and C)
N1. pre-step (non-packed) states: = 65802 (the c2 step count).  95%.
N2. good candidates per state: min = max = 1.  90%.
N3. bad candidates per state: 0 or 1; 1 only when it is the undo; never at k=0.  95%.
N4. states where R is edge-adjacent to a nest rhombus: ~13% (c2 probe, with multiplicity); good-candidate count
    there still exactly 1.  90%.  I do not expect "the good candidate exists only because of a nest rhombus" to be
    meaningful: candidates never contain a nest rhombus, so nest rhombi can only REMOVE candidates.
N5. distinct gash-neighbourhood types (pieces meeting {P,M,Q}, translated, + which nearby unit triangles are
    outside T_n): a few hundred at most.  60%.
B/B: journeys are the same engine, so the same holds (not a B/C fact).
