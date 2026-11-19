"""Inspect closures: which (z,T) closes, where z lands (nest index, cell), sizes of the differences before/after,
and whether the closing third move lands in a nest whose data differs between Mxy and Myx."""
import sys, os, copy
os.chdir('/home/clio/projects/proofs/code-1118c2-Q597')
sys.argv = [sys.argv[0], '60']
exec(open('/home/clio/projects/proofs/code-1118-Q597/widen.py').read().split('st = Counter()')[0])
def full(b): return (tiles_of(b), tuple(frozenset(b.nest[i]['pos']) for i in range(3)))
def mvland(b, z, T, DV):
    bb = copy.deepcopy(b); Y, ij2, s = mig1(bb, 0, z, DV[T]); return bb, Y, ij2
st = Counter(); n = 0
for key, b0 in census():
    for b, DV in states(b0, cap=CAP):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if not x < y: continue
                bxy = mv(mv(b, x, 'B', DV), y, 'B', DV); byx = mv(mv(b, y, 'B', DV), x, 'B', DV)
                if full(bxy) == full(byx): continue
                st[('nests equal before 3rd', full(bxy)[1] == full(byx)[1])] += 1
                st[('tiles differ before 3rd', len(set(tiles_of(bxy)) ^ set(tiles_of(byx))))] += 0
                zs = set(corners(set(bxy.nest[0]['pos']))) & set(corners(set(byx.nest[0]['pos'])))
                for z in zs:
                    for T in 'BC':
                        b1, Y1, c1 = mvland(bxy, z, T, DV); b2, Y2, c2 = mvland(byx, z, T, DV)
                        if full(b1) == full(b2):
                            st[('closing move: dir', T, 'lands nest', Y1, 'cell', 'same' if c1 == c2 else 'diff')] += 1
for k in sorted(st, key=str):
    if st[k]: print(st[k], k)
