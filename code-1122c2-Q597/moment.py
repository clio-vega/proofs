"""(N) via the first moment: per good configuration (target, class, type), the TOTAL readout value change on each side,
sum_k delta(k) (a vector in the other block), from the Theorem AW table.  usage: python3 moment.py"""
import sys, os
from fractions import Fraction as Fr
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
exec(open('/home/clio/projects/proofs/code-1122c2-Q597/sidecfg.py').read().split('st = Counter()')[0])
A_CLASS = frozenset({(0, 1), (0, -1)})   # nest-A class (class 0 = start class): Lambda_0 edge +-(0,1)
for T in ('B', 'C'):
    C = Wd.configs(T)
    for (ck, typ), (R, rem, R2, cre) in sorted(C.items(), key=str):
        d = sdelta(rem, cre); cl = 0 if ck == A_CLASS else 1
        tot = {s: tuple(sum(v[i] for k, v in d.items() if k[0] == s) for i in range(2)) for s in (0, 1)}
        al = tuple(b-a for a, b in zip(cen(R), cen(R2)))
        print(T, 'class', cl, sorted(ck), typ, 'alpha', tuple(str(x) for x in al), 'total side0', tot[0], 'total side1', tot[1])
