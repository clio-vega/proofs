"""Joint local picture: x's steps (other journey Jo, in M) whose tiles meet the region of y's two middles.
Report per config: types and rhombus displacements of those x steps, in units, relative to alpha."""
from replay import *
from collections import Counter, defaultdict
def disp(s): return tuple(str(y-x) for x, y in zip(cen(rh(s[1])), cen(rh(s[2]))))
pat = defaultdict(Counter)
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, A, B, Jo = clean(r); mA, mB, p, q = middles(A, B)
    reg = set().union(*[s[1] | s[2] for s in mA + mB])
    idx = [i for i, s in enumerate(Jo) if (s[1] | s[2]) & reg]
    w = ''.join(s[0] for s in mA)
    xs = ' '.join(Jo[i][0] + str(disp(Jo[i])).replace("'", '') for i in range(idx[0], idx[-1]+1))
    pat[(side, w)][(xs, tuple(idx) == tuple(range(idx[0], idx[-1]+1)))] += 1
for k in sorted(pat):
    print(k)
    for xs, c in pat[k].most_common(4): print('   ', c, xs)
