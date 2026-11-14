"""Shape-level confluence: is the DESTINATION NEST SHAPE after removing alpha\\I a function of I alone (over all removal orders)?
Also count the information: cells where the mosaic is order-dependent but the shape is not."""
from common import *
from collections import Counter, defaultdict
st = Counter()
def explore(b0, dvec, want):
    shapes = defaultdict(set); tiles = defaultdict(set)
    stack = [b0]; seen = set()
    while stack:
        b = stack.pop(); I = frozenset(b.nest[0]['pos']); tk = tuple(tiles_of(b))
        if (I, tk) in seen: continue
        seen.add((I, tk)); tiles[I].add(tk); shapes[I].add(frozenset(b.nest[want]['pos']))
        for x in corners(I):
            bx = copy.deepcopy(b); Y, ij, s = mig1(bx, 0, x, dvec)
            assert Y == want; stack.append(bx)
    return shapes, tiles
for lab, b0 in census():
    A = b0.nest[0]; uA, wA = A['u'], A['w']
    for name, dvec, want in (('B', wA, 1), ('C', uA, 2)):
        sh, ti = explore(copy.deepcopy(b0), dvec, want)
        st[name+' ideals'] += len(sh)
        st[name+' ideals w/ >1 mosaic'] += sum(len(v) > 1 for v in ti.values())
        st[name+' ideals w/ >1 dest shape'] += sum(len(v) > 1 for v in sh.values())
        st[name+' cells shape-confluent'] += all(len(v) == 1 for v in sh.values())
    st['cells'] += 1
print(dict(st))
