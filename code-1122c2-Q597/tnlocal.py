"""(Tn) locally?  After J_o's class-0 b-step (table config), at the image rhombus: for each class-1 configuration, do its tiles
(a) include a tile removed by the b-step (then absent: excluded), (b) overlap-conflict with a created tile (not equal but
same kind & sharing >= 3 vertices for tri / all 4 for sq would be equal; here we test EQUALITY with created = present)?
usage: python3 tnlocal.py"""
import sys, os
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
exec(open('/home/clio/projects/proofs/code-1122c2-Q597/sidecfg.py').read().split('st = Counter()')[0])
A_CLASS = frozenset({(0, 1), (0, -1)})
def base(R): return min(R[1])
STRAIGHT = {'B': 'a1', 'C': 'a0'}
for Tp in ('C', 'B'):
    C = Wd.configs(Tp)
    (R, rem, R2, cre) = next(v for (ck, typ), v in C.items() if ck == A_CLASS and typ == 'b')
    for (ck, typ), (Rc, remc, R2c, crec) in sorted(C.items(), key=str):
        if ck == A_CLASS: continue
        v = tuple(a - b for a, b in zip(base(R2), base(Rc)))
        need = {shift(t, v) for t in remc if t[0] != 'rh'}
        print('Jo ->', Tp, 'next', typ, '(straight)' if typ == STRAIGHT[Tp] else '', '| tiles removed by the b-step:', len(need & set(rem)), '| tiles created by the b-step:', len(need & set(cre)), '| unknown:', len(need - set(rem) - set(cre)))
