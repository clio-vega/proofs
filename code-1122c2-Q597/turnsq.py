"""(T) on the table: J_o toward T' takes its class-0 b-step (config at canonical R, base 0); then at the image rhombus R2
each class-1 T'-configuration (a0, a1, b) translated to R2.  Is (b-square - alpha) created by that next step?
alpha = delta_T, T the clean target.  usage: python3 turnsq.py"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
exec(open('/home/clio/projects/proofs/code-1122c2-Q597/sidecfg.py').read().split('st = Counter()')[0])
A_CLASS = frozenset({(0, 1), (0, -1)})
DELTA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
def base(R): return min(R[1])
for T, Tp in (('B', 'C'), ('C', 'B')):
    ma = tuple(-x for x in DELTA[T]); C = Wd.configs(Tp)
    (R, rem, R2, cre) = next(v for (ck, typ), v in C.items() if ck == A_CLASS and typ == 'b')
    sq = next(u for u in rem if u[0] == 'sq'); target = shift(sq, ma)
    print(T, Tp, 'turn square - alpha created by the b-step itself:', target in cre)
    for (ck, typ), (Rc, remc, R2c, crec) in sorted(C.items(), key=str):
        if ck == A_CLASS: continue
        v = tuple(a - b for a, b in zip(base(R2), base(Rc)))
        if shift(Rc, v) != R2: print('  class-1', typ, 'does not sit at the image rhombus'); continue
        crec_t = {shift(t, v) for t in crec}; remc_t = {shift(t, v) for t in remc}
        print('  next step class 1', typ, ': square-alpha created', target in crec_t, '| removed', target in remc_t)
