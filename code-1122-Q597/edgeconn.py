"""PROVE 11-22: finite check for Lemma U.  In each of the 12 good configurations (both T, both classes, three types),
the four tiles (R + square + 2 triangles) are connected under 'shares a FULL edge' (two common vertices that are
adjacent on both polygons) -- the adjacency the engine's candidate growth uses (Board.neighbours) -- with the rhombus
as root and no step through another rhombus.  Control: the same test after translating one tile by a unit vector
must disconnect it.   usage: python3 edgeconn.py"""
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
def cyc(t):
    hull = poly(t); pts = list(t[1]); out = []
    for z in hull:
        out.append(min(pts, key=lambda P: abs(LE.phi(P) - z)))
    return out
def fedges(t):
    c = cyc(t); return {frozenset((c[i], c[(i+1) % len(c)])) for i in range(len(c))}
def connected(tiles):
    R = next(t for t in tiles if t[0] == 'rh'); seen = {R}; stack = [R]
    while stack:
        a = stack.pop()
        for b in tiles:
            if b not in seen and fedges(a) & fedges(b): seen.add(b); stack.append(b)
    return len(seen) == len(tiles)
st = Counter()
for T in 'BC':
    for c in (0, 1):
        for typ in TYPES:
            R = canon(T, c); rem = step(T, R, typ)[0]
            st[(T, c, typ, 'connected under full-edge adjacency', connected(rem))] += 1
            t = nonrh(rem)[0]; moved = [x if x != t else tr(x, (1, 0, 0, 0)) for x in rem]
            st[('CONTROL one tile shifted: still connected (expect False)', connected(moved))] += 1
for k in sorted(st, key=str): print(st[k], k)
