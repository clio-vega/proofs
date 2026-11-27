# PREDICTION — PROVE 2026-11-27 c1 (written before any run)
Reformulation: Prop wb(2) does not need (Q2). It needs only that G(E_{k+3}) (M0 tiles, removed by J_E) and
G(E^L_{k+3}) lie on ONE tiling. Claim (Q2'): every tile of G(E_{k+3}) is still on E^L's board before E^L_{k+3}.
Then the full-set coupling table of wb.py (bad pairs overlap, census pairs don't) fixes the pair.
(Q2') reduces, by Prop rerun(1)'s template (dinv + prefix minimality + shadow + bsall + fin), to a finite window check:
 W1: G(E_{k+3}) misses G(L_0), G(L_1), G(E^L_j) for 2<=j<=k+2 (local words), all k. PREDICT: holds, window closes at small k.
 S1: R^E_{k+3} (created by E_{k+2} = T_E's (0,b)) is NOT of T_E's packed shape. PREDICT: 60% it is not packed shape
     (then "J_E steps there" is a shape fact). If it IS packed shape, the no-stop needs the nest-ideal condition: harder.
B/B control: the step-0 (a)/(c) censuses; prediction: (Q2) holds at k=0 too, 4 classes same.
# B/B control (written before the run): on B/B pairs passing the same role filter (L at (1,a_T), L^E starts (1,b)),
# (LEN) |mid_E| = k+4 FAILS on (nearly) every record; the stale-set "read later" may be nonzero.
