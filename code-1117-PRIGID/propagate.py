"""Does unit propagation reconstruct S from its pure cells?  Start: pure cells of S placed.  Repeat: find a needed half-edge
with exactly ONE admissible parallelogram (no overlap with placed tiles, no half-edge reuse); place it.  Report completions/stalls."""
import sys, pickle, random
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E
def setup(a, b):
    tiles = E.make_tiles(a, b); poly = [E.ccw(vs) for _, vs in tiles]; idx = {vs: i for i, (_, vs) in enumerate(tiles)}
    allpts = set(v for _, vs in tiles for v in vs); HE = {}
    for i, P in enumerate(poly):
        for j in range(len(P)): HE.setdefault((P[j], P[(j+1) % len(P)]), []).append(i)
    bd = set(h for h in HE if all(E.cross(h[0], h[1], x) >= 0 for x in allpts))
    n = len(tiles); conf = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if E.overlap(poly[i], poly[j]): conf[i] |= 1 << j; conf[j] |= 1 << i
    return tiles, poly, idx, HE, bd, conf
def place(i, poly, needed, used):
    P = poly[i]
    for j in range(len(P)):
        h = (P[j], P[(j+1) % len(P)])
        if h in used: return False
    for j in range(len(P)):
        h = (P[j], P[(j+1) % len(P)]); used.add(h)
        if h in needed: needed.discard(h)
        else: needed.add((h[1], h[0]))
    return True
def propagate(T, S):
    tiles, poly, idx, HE, bd, conf = S
    needed, used, forbid, placed = set(bd), set(), 0, set()
    for vs in T:
        if vs[0] != 'par':
            i = idx[vs[1]]; assert place(i, poly, needed, used); forbid |= conf[i] | (1 << i); placed.add(i)
    while needed:
        forced = None
        for h in needed:
            c = [i for i in HE.get(h, []) if not (forbid >> i) & 1 and all((poly[i][j], poly[i][(j+1) % len(poly[i])]) not in used for j in range(len(poly[i])))]
            if len(c) == 1: forced = c[0]; break
            if len(c) == 0: return 'dead'
        if forced is None: return 'stall'
        place(forced, poly, needed, used); forbid |= conf[forced] | (1 << forced); placed.add(forced)
    return 'ok' if placed == set(idx[vs] for _, vs in T) else 'wrong'
for fn in sys.argv[1:]:
    a, b, M = pickle.load(open(fn, 'rb')); S = setup(a, b); M = list(M); random.seed(1)
    sample = M if len(M) <= 3000 else random.sample(M, 3000)
    print(fn, len(sample), Counter(propagate(T, S) for T in sample))
