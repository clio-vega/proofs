"""Phase 1 step 2: is the diagonal dressing well defined -- before or after?"""
import sys
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine
from bead import maya, from_maya, trim, bilinear

n_moves = 0; disagree = 0
for e in (2, 3, 4, 5, 6):
    for n in range(0, 10):
        for lam in (engine.parts_of(n) if n else ((),)):
            lo = -(len(lam) + 2 * e + 6); M = maya(lam, lo)
            hi = (lam[0] if lam else 0) + 2 * e + 6
            for b in range(lo, hi + 1):
                if b not in M or b + e in M: continue
                r = bilinear(b + e, b, M)
                if r is None: continue
                _, Mp = r
                before = sum(1 for j in range(b + 1, b + e) if j in M)
                after  = sum(1 for j in range(b + 1, b + e) if j in Mp)
                n_moves += 1
                if before != after: disagree += 1
print(f"bead moves examined: {n_moves};  N(before) != N(after): {disagree}")
