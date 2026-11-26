"""Triangle moves of the four journeys: is each triangle moved at most once per journey? Is the displacement a function of (target, side, step type)?
Then compare the symmetric differences  moved(E^L) - moved(E)  vs  moved(L^E) - moved(L)."""
import pickle, glob
from collections import Counter, defaultdict
from lookstep import delta
import sys
sys.argv=['x','5','1','24']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
st = Counter(); dispmap = defaultdict(Counter)
def moves(J):
    M = Counter(); D = defaultdict(lambda: (0, 0)); ty = typed(J)
    for s, t in zip(J, ty):
        for k, v in delta(s).items():
            M[k] += 1; D[k] = tuple(a+b for a, b in zip(D[k], v))
    return M, D
for f in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(f, 'rb')):
        for (T, A, Ar, p, la, lb) in info:
            for J in (A, Ar):
                M, D = moves(J)
                st[('max times a triangle is moved in one journey', max(M.values()))] += 1
                for s, t in zip(J, typed(J)):
                    for k, v in delta(s).items(): dispmap[(T, k[0], s[0])][v] += 1
for k in sorted(st, key=str): print(st[k], k)
for k in sorted(dispmap): print('target,side,steptype', k, dict(dispmap[k]))
