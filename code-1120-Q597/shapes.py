import sys, pickle, glob
from collections import Counter
sys.argv = [sys.argv[0], '5']; exec(open('potential.py').read().split('st = Counter()')[0])
def rhom(S): return next(t for t in S if t[0] == 'rh')
def dirn(vs, i):
    pts = sorted(proj(vs, 1-i)) if False else sorted({(v[:2] if i == 0 else v[2:]) for v in vs})
    a, b = pts; d = (b[0]-a[0], b[1]-a[1]); return max(d, (-d[0], -d[1]))
st = Counter()
for fn in sorted(glob.glob('dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for TT, M0, J in D['BC']:
        for k, (typ, g, n) in enumerate(J):
            R = rhom(g); st[(TT, 's0', dirn(R[1], 0), 's1', dirn(R[1], 1))] += 1
            for t in g:
                if t[0] == 'sq': st[(TT, 'SQ s0', dirn(t[1], 0), 's1', dirn(t[1], 1))] += 1
        st[(TT, 'first step type', J[0][0])] += 1
for k in sorted(st, key=str): print(st[k], k)
