"""PROVE 11-21 c2, Lemma shared-sheet (finite check on the configuration table of Lemma word, code-1120c2-Q597/words.py):
for the class u = +-(0,1) (allowed for both B and C), compare the unique good configuration of each type for B and C:
same removed tiles (hexagon)?  same image rhombus?  same created tiles?
Control: the other classes (+-(1,-1) for B, +-(1,0) for C) have disjoint keys, so no comparison is possible there;
lesion: compare B's class-0 a0 with C's class-0 a1 (must differ).   usage: python3 class0same.py"""
import os, sys
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '.')
exec(open('words.py').read().split('\ndef run(T, L, st):')[0])
CB, CC = configs('B'), configs('C')
k0 = frozenset({(0, 1), (0, -1)})
print('B keys', sorted((sorted(k), t) for k, t in CB)); print('C keys', sorted((sorted(k), t) for k, t in CC))
for typ in ('a0', 'a1', 'b'):
    vb, vc = CB[(k0, typ)], CC[(k0, typ)]
    print('class0', typ, 'R same', vb[0] == vc[0], 'removed same', vb[1] == vc[1], "R' same", vb[2] == vc[2], 'created same', vb[3] == vc[3])
print('LESION B a0 vs C a1 removed same:', CB[(k0, 'a0')][1] == CC[(k0, 'a1')][1])
print('LESION B a1 vs C a0 removed same:', CB[(k0, 'a1')][1] == CC[(k0, 'a0')][1])
print('LESION B b vs C a0 removed same:', CB[(k0, 'b')][1] == CC[(k0, 'a0')][1])
