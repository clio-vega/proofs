"""Q590 well-posed form: is single-rhombus migration out of A CONFLUENT (a growth over the lattice J(alpha))?
For each census mosaic and each target direction d in {N_A (->B), E_A (->C)}: explore ALL outer-corner removal sequences;
state(I) must be independent of the sequence (compare full tilings).  Also: does every single migration land in the
intended nest?  Prints counts of diamonds tested / failed, and whether (ii)'s / (iii)'s orders are among the sequences."""
from common import *
from collections import Counter
st = Counter(); ex = []
def explore(b0, dvec, want):
    memo = {}  # frozenset(ideal) -> tiles
    stack = [b0]; nd = 0; bad = 0; wrong = 0; diam = Counter()
    seen = set()
    while stack:
        b = stack.pop(); I = frozenset(b.nest[0]['pos'])
        key = (I, tuple(tiles_of(b)))
        if key in seen: continue
        seen.add(key)
        if I in memo and memo[I] != key[1]: bad += 1; continue
        memo[I] = key[1]
        for x in corners(I):
            bx = copy.deepcopy(b)
            try: Y, ij, s = mig1(bx, 0, x, dvec)
            except Exception as e: st['ERR '+str(e)[:30]] += 1; continue
            if Y != want: wrong += 1; continue
            stack.append(bx)
    return memo, bad, wrong
n = 0
for lab, b0 in census():
    n += 1; A = b0.nest[0]; uA, wA = A['u'], A['w']
    for name, dvec, want in (('->B (N_A)', wA, 1), ('->C (E_A)', uA, 2)):
        memo, bad, wrong = explore(copy.deepcopy(b0), dvec, want)
        st[name+' ideals'] += len(memo); st[name+' conflicts'] += bad; st[name+' wrong-nest'] += wrong
        st[name+' cells confluent'] += (bad == 0 and wrong == 0)
        if bad and len(ex) < 3: ex.append((lab, name))
    st['cells'] = n
    if n % 20 == 0: print(n, dict(st), flush=True)
print(dict(st)); print(ex)
