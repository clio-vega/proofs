"""FIRST MOVE, as briefed: the nest {Sigma^{<=j}} for witness A and witness B
   of thm:blind.  PLUS the prior question the brief did not ask: are the two
   witnesses realisable AT ALL in the class (H3) carves out?

   The nest lives in y-space.  A certificate is realisable in the abstract class
   iff there are g>=0, a box B and sigma with { w(y) : y in B cap {sum y = sigma} }
   equal to the certificate's width multiset (all widths >= 1).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ycoord import wmap, defect, centred_sum, pos
from gen import is_pf2
from itertools import product
from collections import Counter

A = [(1,1),(2,2),(1,5),(2,6)]   # sum (1,3,4,6,4,3,1)  NOT PF2
B = [(1,1),(2,2),(3,3),(4,4)]   # sum (1,3,6,10,6,3,1)  PF2

def report(name, ws):
    m = len(ws[0])
    half = [(sum(w)-m)//2 for w in ws]
    c = max(half)
    k = [c - h for h in half]            # k(nu) = c - Lambda(nu)
    print(f"  {name}: widths {ws}")
    print(f"       half-widths {half}, defect levels k = {k}")
    for j in range(0, c+1):
        lev = [w for w, kk in zip(ws, k) if kk <= j]
        print(f"       |Sigma^<=({j})| = {len(lev)}   widths {sorted(lev)}")
    print(f"       sum {centred_sum(ws)}  PF2={is_pf2(centred_sum(ws))}")

print("=== nest data for the two thm:blind witnesses (from the width vectors) ===")
report("A (not PF2)", A)
report("B (PF2)", B)

print()
print("=== PRIOR QUESTION: is either witness realisable in the abstract class? ===")
print("m=2 exhaustive search over g_1,g_2 <= 12, boxes in [-12,12], all sigma.")
targetA, targetB = sorted(A), sorted(B)
hitsA = hitsB = 0
R = range(-12, 13)
for g1 in range(0, 13):
    for g2 in range(0, 13):
        g = (g1, g2)
        for P1 in R:
            for Q1 in range(P1, 13):
                for P2 in R:
                    for Q2 in range(P2, 13):
                        for sig in range(P1+P2, Q1+Q2+1):
                            Y = [(y1, sig-y1) for y1 in range(P1, Q1+1)
                                 if P2 <= sig-y1 <= Q2]
                            if len(Y) != 4: continue
                            ws = [wmap(y, g) for y in Y]
                            if any(min(w) <= 0 for w in ws): continue
                            s = sorted(ws)
                            if s == targetA: hitsA += 1
                            if s == targetB: hitsB += 1
print(f"  realisations of witness A: {hitsA}")
print(f"  realisations of witness B: {hitsB}")

print()
print("=== why: injectivity of y -> w_1-w_2 at m=2 ===")
print("  w_1-w_2 = (g_1-g_2) - y_1 + y_2  and  y_1+y_2 = sigma, so y |-> w_1-w_2")
print("  is injective on a slice.  Width DIFFERENCES of the certificates:")
for nm, ws in (("A", A), ("B", B)):
    print(f"    {nm}: {[w[0]-w[1] for w in ws]}  -> distinct? {len(set(w[0]-w[1] for w in ws))==len(ws)}")
