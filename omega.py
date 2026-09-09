"""
ENGINE 2: direct Young-diagram border strips (code-disjoint from the abacus).
Test the candidate operator identity   omega R_e(t) omega = t^(e-1) R_e(1/t).
"""
import sys
sys.path.insert(0, '/home/clio/projects/reviews/2026-09-08-selfreview-code')
from ribbon import border_strips, parts, cells      # cross-checked 764/764 vs Q105 Thm A

def conj(lam):
    if not lam: return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))

# ---- Step 1: does transposition complement the height of a border strip? ----
bad = 0; tot = 0; hist = {}
for e in (1,2,3,4,5,6):
    for lam in parts(0):
        pass
    for N in range(0, 10):
        for lam in parts(N):
            for (mu, h) in border_strips(lam, e):     # mu/lam a strip of size e
                lc, mc = conj(lam), conj(mu)
                # find the transposed strip's height
                strips = dict(border_strips(lc, e))
                tot += 1
                if mc not in strips:
                    bad += 1; print("NOT A STRIP:", lam, mu, e); continue
                h2 = strips[mc]
                hist.setdefault(e, set()).add((h, h2))
                if h + h2 != e - 1:
                    bad += 1
                    print(f"HEIGHT MISMATCH e={e} {lam}->{mu}: h={h} h'={h2}")
print(f"transpose test: {tot} (lam, strip) pairs, {bad} failures")
for e in sorted(hist):
    print(f"   e={e}: observed (h, h') pairs = {sorted(hist[e])}   (need h+h'={e-1})")
