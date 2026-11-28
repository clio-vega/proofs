"""FC2: two rhombi in flight.  For a state b and corners x != y of nest A (any pair, incl. the SW pair after x leaves):
x takes j steps toward T, then y takes k steps toward T' (on that state), then FINISH x then y  vs  FINISH y then x.
Report agreement over all (j,k) where both finishes are defined.  T/T' = B/C (claim) and B/B (control).
usage: python3 fc2.py NSEED JMAX KMAX"""
import sys
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
def packed(b, R, src=0):
    vs = b.tiles[R][1]
    for X, info in enumerate(b.nest):
        if X == src: continue
        ij = b.packed_at(info, vs)
        if ij is not None and ij not in info['pos'] and b.is_ideal(set(info['pos']) | {ij}): return True
    return False
def steps(b, R, d, n):
    for _ in range(n):
        if packed(b, R): return None
        R = b.step(R, d)
        if R is None: return None
    return R
def finish(b, R, d):
    b.journey(R, 0, d)
def trial(b, x, y, dx, dy, j, k, adj):
    b = copy.deepcopy(b); d1, d2 = unit_dir(dx), unit_dir(dy)
    Rx = b.nest[0]['pos'].pop(x)
    Rx = steps(b, Rx, d1, j)
    if Rx is None: return 'x-pre'
    Ry = b.nest[0]['pos'].pop(y)
    Ry = steps(b, Ry, d2, k)
    if Ry is None: return 'y-pre'
    out = []
    for order in (0, 1):
        c = copy.deepcopy(b)
        try:
            if order == 0: finish(c, Rx, d1); finish(c, Ry, d2)
            else: finish(c, Ry, d2); finish(c, Rx, d1)
        except RuntimeError: out.append(None); continue
        out.append(tuple(tiles_of(c)))
    if None in out: return 'finish-err'
    return out
st = Counter(); NSEED, JM, KM = map(int, sys.argv[1:4])
for kk, (key, b0) in enumerate(census()):
    if kk >= NSEED: break
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        pairs = [(x, y, 'MC') for x in cs for y in cs if x != y]
        pairs += [(x, (x[0], x[1]-1), 'SW') for x in cs if (x[0], x[1]-1) in corners(S - {x})]
        for x, y, kind in pairs:
            for T, T2, lab in (('B', 'C', 'BC'), ('C', 'B', 'CB'), ('B', 'B', 'BB')):
                if kind == 'SW' and lab == 'BB': continue
                ref = None
                for j in range(JM + 1):
                    for k in range(KM + 1):
                        r = trial(b, x, y, dirs[T], dirs[T2], j, k, kind)
                        if isinstance(r, str): st[(kind, lab, 'undef', r)] += 1; continue
                        st[(kind, lab, 'defined')] += 1
                        st[(kind, lab, 'orders agree')] += r[0] == r[1]
for k in sorted(st, key=str): print(st[k], k)
