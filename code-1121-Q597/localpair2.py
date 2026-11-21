"""PROVE 11-21: the square-sharing table with offsets and rhombus coincidences.  For every (c_J, c_y, v) with a common
SQUARE (c_y translated by v): y's rhombus R_y+v equal to Jo's R_J (before) or R'_J (after the step)?  Also: do the two
removed sets share a triangle too, and is the y hexagon equal to the Jo hexagon?"""
import sys, os
sys.argv = ['x', '0', '0']
exec(open('/home/clio/projects/proofs/code-1121-Q597/twoword.py').read().split('st = Counter()')[0])
from collections import Counter
NA = frozenset({(0, 1), (0, -1)})
def kidx(ck): return 0 if ck == NA else 1
for Ty, TJ in (('B', 'C'), ('C', 'B')):
    for (T1, ckJ, tJ), (RJ, remJ, RnJ, creJ) in sorted(CF.items(), key=str):
        if T1 != TJ: continue
        for (T2, cky, ty), (Ry, remy, Rny, crey) in sorted(CF.items(), key=str):
            if T2 != Ty: continue
            for u in remJ:
                if u[0] != 'sq': continue
                for w in remy:
                    if w[0] != 'sq': continue
                    v = tuple(a - b for a, b in zip(min(u[1]), min(w[1])))
                    if tr(w, v) != u: continue
                    Ryv = tr(Ry, v); remyv = set(tr(t, v) for t in remy)
                    tri = sorted(t[0] for t in remyv & set(remJ) if t[0] == 'tri')
                    print(Ty, 'Jo', (kidx(ckJ), tJ), 'y', (kidx(cky), ty), 'v', v, 'R_y=R_J', Ryv == RJ, "R_y=R'_J", Ryv == RnJ,
                          'shared tris', len(tri), 'same hexagon', remyv == set(remJ) and Ryv == RJ)
