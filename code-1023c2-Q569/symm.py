"""Action of the dihedral group D_6 of the triangular lattice on the tiles.

CCW traversal of sigma(P) for an orientation-REVERSING sigma (d -> m-d):
the image vertex cycle is CW, so reverse it; reversing a closed walk
requires reversing the ORDER *and* negating each direction (d -> d+3).
Edge i keeps its label l_i under both operations.
    out = [((m - d + 3) % 6, l) for (d,l) in reversed(steps)]
(An earlier version reversed the order only.  That merely permutes the edges
into a different polygon -- which still closes and is still CCW, so neither
the closure check nor the area check detects it.)
"""
from puzzles import TILES, D, ccw_area2, closes, verts, up_tri, dn_tri

def canon(steps):
    n=len(steps); best=None
    for i in range(n):
        r=tuple(steps[i:]+steps[:i])
        if best is None or r<best: best=r
    return best

def act(steps, g):
    t,p = g
    if t=='r':
        out=[((d+p)%6,l) for d,l in steps]
    else:
        out=[((p-d+3)%6,l) for d,l in reversed(steps)]
    assert closes(out), (g,steps)
    assert ccw_area2(out)==ccw_area2(steps), ('area changed',g,ccw_area2(steps),ccw_area2(out))
    assert sorted(l for _,l in out)==sorted(l for _,l in steps), 'labels changed'
    return out

GROUP=[('r',j) for j in range(6)]+[('s',m) for m in range(6)]
NAME={('r',j):('id' if j==0 else 'rot%d'%(60*j)) for j in range(6)}
for m in range(6): NAME[('s',m)]='refl_%d'%m
INVG={('r',j):('r',(-j)%6) for j in range(6)}
for m in range(6): INVG[('s',m)]=('s',m)

CAN={nm:canon(st) for nm,st in TILES.items()}
INV={v:k for k,v in CAN.items()}
assert len(INV)==len(CAN), "two tiles share a canonical form!"
def image(nm,g): return INV.get(canon(act(TILES[nm],g)))

KTW=['zero_A','zero_B','one_A','one_B','rh_014','rh_230','rh_125']
