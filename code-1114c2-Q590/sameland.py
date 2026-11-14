"""From each reachable state (all removal orders), do ALL outer corners of the current A-ideal land at the SAME destination box?"""
from common import *
from collections import Counter
st = Counter()
for key, b0 in census():
    A = b0.nest[0]; uA, wA = A['u'], A['w']
    for name, dvec, want in (('B', wA, 1), ('C', uA, 2)):
        stack = [b0]; seen = set(); okcell = True
        while stack:
            b = stack.pop(); I = frozenset(b.nest[0]['pos']); tk = tuple(tiles_of(b))
            if (I, tk) in seen: continue
            seen.add((I, tk)); cs = corners(I); lands = set()
            for x in cs:
                bx = copy.deepcopy(b); Y, ij, s = mig1(bx, 0, x, dvec); lands.add(ij); stack.append(bx)
            if len(cs) >= 2:
                st[name+' states with >=2 corners'] += 1
                st[name+' ... all corners same landing'] += len(lands) == 1
                okcell &= len(lands) == 1
        st[name+' cells where it always holds'] += okcell
    st['cells'] += 1
print(dict(st))
