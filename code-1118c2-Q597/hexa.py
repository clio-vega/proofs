"""Retreat test (a), DREAM 11-18: for each NON-commuting B/B pair (x,y) in the <=5 census (cap 60 = full census),
does a third move close the square into a hexagon?  The two final states can agree only if nest A loses the same cells,
so z = z' is forced.  Test: exists corner z of A minus {x,y} (a corner in BOTH Mxy and Myx) and target T in {B,C}
with Mxy.z == Myx.z (tiles and nests).  Also report, for each z,T, whether the two journeys of z re-synchronise."""
import sys, os, copy
os.chdir('/home/clio/projects/proofs/code-1118c2-Q597')
sys.argv = [sys.argv[0], '60']
exec(open('/home/clio/projects/proofs/code-1118-Q597/widen.py').read().split('st = Counter()')[0])
def full(b): return (tiles_of(b), tuple(frozenset(b.nest[i]['pos']) for i in range(3)))
st = Counter(); n = 0
for key, b0 in census():
    for b, DV in states(b0, cap=CAP):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if not x < y: continue
                bxy = mv(mv(b, x, 'B', DV), y, 'B', DV); byx = mv(mv(b, y, 'B', DV), x, 'B', DV)
                if full(bxy) == full(byx): continue
                n += 1
                zs = set(corners(set(bxy.nest[0]['pos']))) & set(corners(set(byx.nest[0]['pos'])))
                closes = []
                for z in zs:
                    for T in 'BC':
                        try:
                            if full(mv(bxy, z, T, DV)) == full(mv(byx, z, T, DV)): closes.append((z, T))
                        except Exception as e:
                            st['error ' + type(e).__name__] += 1
                st[('closes after 3rd move', bool(closes))] += 1
                st[('#closing (z,T)', len(closes), 'of', 2*len(zs))] += 1
print('noncommuting B/B pairs', n)
for k in sorted(st, key=str): print(st[k], k)
