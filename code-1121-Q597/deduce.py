"""PROVE 11-21: mechanical check of the table step of Prop. CS.  From localpair2.out (squares) and localpair3-tri.out
(triangles): Y_sq = y-configs sharing a square with a Jo step in S_sq, Y_tri = sharing a triangle with a Jo step in S_tri,
rows with R_y = R_J or R_y = R'_J removed (excluded by Lemma P2).  Hypotheses (weakest forms used in the tex):
  T_y=B: S_sq = {(0,a0),(1,b)}, S_tri = {(0,a0),(0,a1)}  (Jo->C);   T_y=C: S_sq = {(0,a1),(1,b)}, S_tri = {(0,a0),(0,a1)}.
Assert Y_sq ∩ Y_tri = {straight}.  Lesion: S_sq with (1,a0)/(1,a1) added (the model's violating Jo step) must break it."""
import re
def rows(fn):
    out = []
    for line in open(fn):
        m = re.match(r"(\w) Jo \((\d), '(\w+)'\) y \((\d), '(\w+)'\) v \(.*?\) R_y=R_J (\w+) R_y=R'_J (\w+)", line)
        if m:
            Ty, kJ, tJ, ky, ty, c1, c2 = m.groups()
            out.append((Ty, (int(kJ), tJ), (int(ky), ty), c1 == 'True' or c2 == 'True'))
    return out
SQ, TRI = rows('localpair2.out'), rows('localpair3-tri.out')
STRAIGHT = {'B': (1, 'a1'), 'C': (1, 'a0')}
HYP = {'B': ({(0, 'a0'), (1, 'b')}, {(0, 'a0'), (0, 'a1')}), 'C': ({(0, 'a1'), (1, 'b')}, {(0, 'a0'), (0, 'a1')})}
LES = {'B': (1, 'a0'), 'C': (1, 'a1')}
for Ty in 'BC':
    for lesion in (False, True):
        Ssq, Stri = set(HYP[Ty][0]), set(HYP[Ty][1])
        if lesion: Ssq.add(LES[Ty])
        Ysq = {y for (T, J, y, co) in SQ if T == Ty and J in Ssq and not co}
        Ytri = {y for (T, J, y, co) in TRI if T == Ty and J in Stri and not co}
        res = Ysq & Ytri
        print(Ty, 'lesion' if lesion else 'hyp', 'Y_sq', sorted(Ysq), 'Y_tri', sorted(Ytri), '->', sorted(res),
              'OK' if res == {STRAIGHT[Ty]} else 'NOT ONLY STRAIGHT')
print('rows', len(SQ), len(TRI))
