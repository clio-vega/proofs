"""PROVE 11-21: local square-sharing table (board-free, exhaustive).  For each pair (T_y, T_J) in {(B,C),(C,B)}, every good
configuration c_J of T_J (class k_J rel. the nest-A class, type t_J), and every good configuration c_y of T_y placed at
ANY translate such that the two removed sets share a SQUARE: which (k_y, t_y) occur?  Optional consistency: the two
hexagons' other non-rhombus tiles must be equal or interior-disjoint (necessary if both are M0 tiles) [column 'cons'].
Class index: 0 = the nest-A class (0,1)/(0,-1), 1 = the other allowed class."""
import sys, os
sys.argv = ['x', '0', '0']
exec(open('/home/clio/projects/proofs/code-1121-Q597/twoword.py').read().split('st = Counter()')[0])
from collections import Counter, defaultdict
NA = frozenset({(0, 1), (0, -1)})
def kidx(ck): return 0 if ck == NA else 1
st = Counter()
for Ty, TJ in (('B', 'C'), ('C', 'B')):
    for (T1, ckJ, tJ), (RJ, remJ, RnJ, creJ) in CF.items():
        if T1 != TJ: continue
        for (T2, cky, ty), (Ry, remy, Rny, crey) in CF.items():
            if T2 != Ty: continue
            # every translate that makes a square of remy coincide with a square of remJ
            for u in remJ:
                if u[0] != 'sq': continue
                for w in remy:
                    if w[0] != 'sq': continue
                    v = tuple(a - b for a, b in zip(min(u[1]), min(w[1])))
                    if tr(w, v) != u: continue
                    remy_t = [tr(t, v) for t in remy]; Ry_t = tr(Ry, v)
                    others_y = [t for t in remy_t if t != u]; others_J = [t for t in remJ if t != u]
                    cons = all(s == t or not overlap(s, t) for s in others_y + [Ry_t] for t in others_J + [RJ]) 
                    st[(Ty, 'Jo', (kidx(ckJ), tJ), 'y', (kidx(cky), ty), 'cons(all tiles incl. rhombi)', cons)] += 1
for k in sorted(st, key=str): print(st[k], k)
