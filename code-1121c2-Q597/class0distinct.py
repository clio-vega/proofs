"""PROVE 11-21 c2: the three class-0 good configurations (a0, a1, b) have pairwise distinct hexagons, for B and for C
(needed in Lemma lock-step: the hexagon determines the type).   usage: python3 class0distinct.py"""
import os, sys
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '.')
exec(open('words.py').read().split('\ndef run(T, L, st):')[0])
k0 = frozenset({(0, 1), (0, -1)})
for T in 'BC':
    C = configs(T); hs = {t: C[(k0, t)][1] for t in ('a0', 'a1', 'b')}
    print(T, 'class-0 hexagons pairwise distinct:', len(set(hs.values())) == 3)
