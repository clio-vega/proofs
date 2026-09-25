"""w(R) := <s_R, p_r> = the coefficient of p_r/z_r in s_R, r=|R|.
For nu=(r) the only ribbon tableau is the single step (empty -> R), so ANY
Murnaghan-Nakayama-type rule for cylindric Schur functions must have
   weight of the step R  =  w(R).
This is exactly the quantity AKO's prop:stackedRibbon evaluates."""
from cyl import *

def w(R):
    return c_nu(R)[(R.size(),)]

def geom(R):
    """structural data of a cylindric diagram: vertical edges, 2x2 blocks, connectivity, rows."""
    cells = set(R.cells())
    x, y, n = R.x, R.y, R.x + R.y
    vert = 0; sq = 0
    for (i, c) in cells:
        if R.inside(i+1, c): vert += 1
        if R.inside(i, c+1) and R.inside(i+1, c) and R.inside(i+1, c+1): sq += 1
    # connectivity on the cylinder
    adj = {}
    for (i, c) in cells:
        nb = []
        for (j, d) in [(i, c+1), (i, c-1), (i+1, c), (i-1, c)]:
            if R.inside(j, d): nb.append(R.rep(j, d))
        adj[(i, c)] = nb
    seen = set(); stack = [next(iter(cells))]
    while stack:
        v = stack.pop()
        if v in seen: continue
        seen.add(v); stack.extend(adj[v])
    conn = (len(seen) == len(cells))
    return dict(size=len(cells), vert=vert, sq=sq, conn=conn, n=n, x=x, y=y)
