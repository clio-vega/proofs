"""Which readout side does each step type touch?  sdelta(step) (bs6.py): per-key change of the pure-triangle readout.
Also: mixed-difference check  Xy-X == Yx-Y  (as dicts) on every dumped pair.   usage: python3 stepside.py SIZE"""
import sys, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
def sub(D, E):
    out = dict(D)
    for k, v in E.items():
        w = tuple(a-b for a, b in zip(out.get(k, (0, 0)), v)); out[k] = w
        if w == (0, 0): del out[k]
    return out
st = Counter()
for P in PAIRS:
    d1 = sub(P['Xy'], P['X']); d2 = sub(P['Yx'], P['Y'])
    st[('MIXED Xy-X == Yx-Y', P['kind'], d1 == d2)] += 1
    st[('MIXED supp == I', P['kind'], set(d1) == P['I'])] += 1
    Ts = ('B', 'C') if P['kind'] == 'BC' else ('B', 'B')
    for (t, A, Bq, Jo), T in zip(P['roles'], Ts):
        for J in (A, Bq):
            TJ = typed(J)
            for s, ct in zip(J, TJ):
                d = sdelta(s)
                if d is None: st[('STEP', T, ct, 'KEYS NOT PRESERVED')] += 1; continue
                st[('STEP', T, ct, 'sides', tuple(sorted({k[0] for k in d})), 'n', len(d))] += 1
for k in sorted(st, key=str): print(st[k], k)
