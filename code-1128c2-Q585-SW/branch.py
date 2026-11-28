"""Branch geometry of (SW): positions as gash segments; classes; step relations around the branch.  usage: branch.py S0 S1"""
import sys
src = open('/home/clio/projects/proofs/code-1128c2-Q585-SW/gpath.py').read().split("N0, N1 =")[0]
exec(src)
from fractions import Fraction as Fr
def cls(p):
    d = (p[2][0]-p[0][0], p[2][1]-p[0][1]); d = (Fr(d[0], 2), Fr(d[1], 2)); return d if d > (0, 0) else (-d[0], -d[1])
def mid(p): return p[1]
def sub(a, b): return (a[0]-b[0], a[1]-b[1])
st = Counter(); S0, S1 = int(sys.argv[1]), int(sys.argv[2])
for k, (key, b0) in enumerate(census()):
    if k < S0 or k >= S1: continue
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            b1 = copy.deepcopy(b); K, _ = path(b1, x, dirs['B']); G, _ = path(b1, y, dirs['C'])
            b2 = copy.deepcopy(b); H, _ = path(b2, x, dirs['C']); Gp, _ = path(b2, y, dirs['B'])
            t = next(i for i in range(len(K)) if K[i] != H[i]); D = sub(mid(K[0]), mid(G[0]))
            f = G[t]
            st[('class f', cls(f), 'class H[t+1]', cls(H[t+1]))] += 1
            st[('f - shadowend', sub(mid(f), mid(G[t-1])) == sub(mid(K[t-1]), mid(K[t-2])) if t >= 2 else None)] += 1
            st[('R: K[t]-(K-D)[t] == D', sub(mid(Gp[t+1]), mid(Gp[t])) == D)] += 1
            st[('H[t+1] - H[t] == D?', sub(mid(H[t+1]), mid(H[t])) == D, 'H[t+1]-f', sub(mid(H[t+1]), mid(f)) == D)] += 1
for k in sorted(st, key=str): print(st[k], k)
