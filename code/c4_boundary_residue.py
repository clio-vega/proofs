#!/usr/bin/env python3
"""Follow-up to c4_boundary_sweep: the MIN v2 per residue class (the governing
quantity for the closest-to-tie), the sharp b-even N1 bound, and the F_k threshold."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'threerow-boundary'))
from c4_boundary_sweep import M_lam, N_lam, hook_f, v2, box_interior, c
from math import comb
from collections import defaultdict

# ---- min v2(M_{b+2}) per (b%8) and per (a%8,b%8,P%8) ----
minv_b = defaultdict(lambda: 99)
minv_class = defaultdict(lambda: 99)
# sharp N1 bound by b-parity
minN1 = {0: 99, 1: 99}
minN2 = 99; minN3_beven = 99
# global min Delta and shape
gmin = {i: (99, None) for i in (1, 2, 3, 4)}

for m in range(8, 401):
    for b in range(c, m):
        a = 2*m - c - b
        if a < b or not box_interior(a, b):
            continue
        P = m - b - c
        if P < 0:
            continue
        M0 = hook_f((a, b, c)); val0 = 2*v2(M0)
        # N bounds
        n1 = v2(N_lam[1](a, b)); n2 = v2(N_lam[2](a, b)); n3 = v2(N_lam[3](a, b))
        minN1[b % 2] = min(minN1[b % 2], n1)
        minN2 = min(minN2, n2)
        if b % 2 == 0:
            minN3_beven = min(minN3_beven, n3)
        # M_{b+2}
        Mb2 = M_lam[2](a, b)
        if Mb2 != 0:
            vv = v2(Mb2)
            minv_b[b % 8] = min(minv_b[b % 8], vv)
            minv_class[(a % 8, b % 8, P % 8)] = min(minv_class[(a % 8, b % 8, P % 8)], vv)
        # Delta
        for i in (1, 2, 3, 4):
            j = b + i
            if j > m:
                continue
            Mf = M_lam[i](a, b)
            if Mf == 0:
                continue
            d = (j + 2*v2(comb(m, j)*Mf)) - val0
            if d < gmin[i][0]:
                gmin[i] = (d, (a, b, m, P))

print("=== sharp content facts (m<=400) ===")
print(f"  min v2(N1): b-even={minN1[0]}  b-odd={minN1[1]}   (memory claimed >=3 b-even, >=2 b-odd)")
print(f"  min v2(N2): {minN2}   (memory claimed >=2)")
print(f"  min v2(N3) b-even: {minN3_beven}   (memory claimed >=1)")
print()
print("=== min v2(M_{b+2}) by b%8 (smallest = closest to M odd) ===")
for k in range(8):
    print(f"  b%8={k}: min v2 = {minv_b[k]}")
print()
print("=== classes (a%8,b%8,P%8) FORCING v2(M_{b+2}) >= 1 (M always even) ===")
forced = sorted(k for k, v in minv_class.items() if v >= 1)
canzero = sorted(k for k, v in minv_class.items() if v == 0)
print(f"  {len(forced)} forcing-even classes, {len(canzero)} classes where M can be odd")
print(f"  forcing-even (a%8,b%8,P%8): {forced}")
print()
print("=== global min Delta(b+i) and extremal shape ===")
for i in (1, 2, 3, 4):
    print(f"  Delta(b+{i}): min={gmin[i][0]} at (a,b,m,P)={gmin[i][1]}")

# ---- F_k threshold: among Q consecutive integers starting at P+1, can we drop
#      k designated factors and keep the product even?  Find sharp Q for F_3. ----
print()
print("=== F_k threshold (min #factors-of-2 in a run of Q consecutive ints, minus k) ===")
def v2run(start, Q):
    return sum(v2(start + t) for t in range(Q))
for Q in range(3, 12):
    # worst case over starts: choose to drop the k largest-v2 members
    worst = {1: 99, 2: 99, 3: 99}
    for start in range(0, 400):
        vs = sorted((v2(start + t) for t in range(Q)), reverse=True)
        for k in (1, 2, 3):
            rem = sum(vs[k:])  # drop k largest
            worst[k] = min(worst[k], rem)
    print(f"  Q={Q}: after dropping k largest-v2, min total v2 -> "
          f"F1(k=1)={worst[1]} F2(k=2)={worst[2]} F3(k=3)={worst[3]}")
