# PREDICTION — PROVE 2026-11-28 c2, (SW) on the uncapped closure. Written 14:57 UTC before any run.

Capped run (swap.py, cap 60 states/cell): 987 column pairs, 987/987.
Full size-5 closure: 5194 states over 134 seeds (mean ~39/seed), so the cap of 60 WAS probably reached in some cells.
Predictions:
1. Size-5 full closure: column-swap pairs tested ~1300 (range 1000-1800); same tiling 100%; B-box and C-box equal 100%.
2. Row swap (mirror): ~1900 tested, 100%.
3. Merge/rejoin: y's 2nd journey ends with a non-empty suffix of x's SOLO journey to the same nest:
   B-side (y->B after x->C): 100%; C-side (y->C after x->B): >=97% (crossing tex: 6 landed before merging on 40 cells x 30 states).
4. Lesion. The brief's (c-1,r-1) is IMPOSSIBLE: an ideal containing x=(c,r) contains (c-1,r) and (c,r-1), so (c-1,r-1)
   is never a corner of a\x. Used instead: the vertex-adjacent corner y'=(c+1,r-1) (both x,y' corners of a).
   Prediction: x->B;y'->C vs x->C;y'->B gives the same tiling in ~0% (it is a subset of the 0/779 control).
5. Size-6 (census_sized chunks 0,1 of 24; never used for SW, used for MC/NS before): 100% hold.
