"""PROVE 11-21 c2 Step 2: does (Sheet) + (P0) give the straight step (i.e. replace (P1))?  Reuses deduce.py's parser on
../code-1121-Q597/localpair2.out (squares) and localpair3-tri.out (triangles).  (Sheet): Jo's removing step is not a
class-1 a-step, S = {(0,a0),(0,a1),(0,b),(1,b)}.  Also each (sq cfg, tri cfg) pair separately, and the c1 hypotheses
as regression (must print OK).  usage: python3 deduce2.py"""
import os, re, itertools
D = '/home/clio/projects/proofs/code-1121-Q597/'
src = open(D + 'deduce.py').read().split('SQ, TRI = ')[0]
exec(src)
SQ, TRI = rows(D + 'localpair2.out'), rows(D + 'localpair3-tri.out')
STRAIGHT = {'B': (1, 'a1'), 'C': (1, 'a0')}
SHEET = {(0, 'a0'), (0, 'a1'), (0, 'b'), (1, 'b')}
ALL = SHEET | {(1, 'a0'), (1, 'a1')}
def Y(tab, Ty, S): return {y for (T, J, y, co) in tab if T == Ty and J in S and not co}
print('rows', len(SQ), len(TRI))
for Ty in 'BC':
    reg = {'B': ({(0, 'a0'), (1, 'b')}, {(0, 'a0'), (0, 'a1')}), 'C': ({(0, 'a1'), (1, 'b')}, {(0, 'a0'), (0, 'a1')})}[Ty]
    r = Y(SQ, Ty, reg[0]) & Y(TRI, Ty, reg[1])
    print(Ty, 'REGRESSION c1 (P1):', sorted(r), 'OK' if r == {STRAIGHT[Ty]} else 'FAIL')
    r = Y(SQ, Ty, SHEET) & Y(TRI, Ty, SHEET)
    print(Ty, 'SHEET+P0 :', 'Y_sq', sorted(Y(SQ, Ty, SHEET)), 'Y_tri', sorted(Y(TRI, Ty, SHEET)), '->', sorted(r))
    r = Y(SQ, Ty, ALL) & Y(TRI, Ty, ALL)
    print(Ty, 'NO HYP   :', sorted(r))
    for js in sorted(ALL):
        for jt in sorted(ALL):
            r = Y(SQ, Ty, {js}) & Y(TRI, Ty, {jt})
            if r: print('  ', Ty, 'sq by', js, 'tri by', jt, '-> clean', sorted(r), 'SHEET-allowed' if js in SHEET and jt in SHEET else 'sheet-forbidden')
