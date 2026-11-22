"""Is the re-run middle a permutation of the middle (identical steps, as (type, removed, created))?  Per role, per side.
usage: python3 perm.py SIZE"""
import sys, pickle
from collections import Counter
SIZE=int(sys.argv[1]); sys.argv=[sys.argv[0],str(SIZE),'1','1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = {'BC': ('a1', 'a0'), 'BB': ('a1', 'a1')}
st = Counter()
for P in PAIRS:
    kind = P['kind']
    sides = {k[0] for k in P['I']}
    side = 'BB' if kind == 'BB' else ('empty' if not sides else int(next(iter(sides))))
    info = []
    for i, (t, A, Bq, Jo) in enumerate(P['roles']):
        mA, mB, p, q = middles(A, Bq)
        info.append((typed(A)[p], typed(Bq)[p], Counter(mA) == Counter(mB), len(mA), len(mB), bs_test(mA, mB),
                     ws(typed(A)[p:len(A)-q]), ws(typed(Bq)[p:len(Bq)-q])))
    L = [i for i in (0, 1) if info[i][0] == (1, aT[kind][i])][0]
    for i in (0, 1):
        rl = 'L' if i == L else 'E'
        st[(kind, side, rl, 'rerun', info[i][1][1], 'perm', info[i][2], 'BS', info[i][5], 'len eq', info[i][3] == info[i][4])] += 1
        if side == 'empty' and rl == 'L': st[('EMPTY L middles', info[i][6], '->', info[i][7])] += 1
for k in sorted(st, key=str): print(st[k], k)
