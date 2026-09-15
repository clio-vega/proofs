import sys; sys.path.insert(0,'/home/clio/projects/proofs/code/2026-09-15-c1-Q149')
from controls import *
print('=== CONTROL 2: realisability of every occupancy word ===')
realisability(emax=6, nmax=14)
print('\n=== CONTROL 4: e=f must yield zero equations ===')
same_rank([1,2,3,4])
print('\n=== CONTROL 1: calibration to Q147 height-only locus ===')
calibration([(1,3),(2,3),(1,4),(2,4),(2,5),(3,4)])
print('\n=== CLASSIFICATION: exhaustive over F_p (includes degenerate strata) ===')
for args in [(1,3,5),(2,3,5),(1,4,3),(2,4,3),(1,5,3)]:
    exhaustive(*args); sys.stdout.flush()
