"""The other journey (q-side: x in My; p-side: y in Mx). Match each middle step of the second journey to a translate of
some M-middle step, forward (same typ, same gone->new) or BACKWARD (gone/new swapped). Report the matching pattern."""
from collections import Counter
from q597lib import *
def shp(step):
    o = cen(rh(step[1])); return tr(step, tuple(-x for x in o))
def rev(step): return (step[0], step[2], step[1])
def shpr(step):
    s = rev(step); o = cen(rh(s[1])); return tr(s, tuple(-x for x in o))
st = Counter(); pats = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side = r['Us'][0][0]
    A, B = (r['Jx'], r['Jxy']) if side == 0 else (r['Jy'], r['Jyx'])
    mA, mB, p, q = middles(A, B)
    fw = {shp(s): i for i, s in enumerate(mA)}; bw = {shpr(s): i for i, s in enumerate(mA)}
    pat = tuple(('F%d' % fw[shp(s)]) if shp(s) in fw else (('R%d' % bw[shp(s)]) if shp(s) in bw else '?') for s in mB)
    pats[(side, r['Us'][0][1], len(mA), pat)] += 1
    st[(side, 'all matched', '?' not in pat)] += 1
for k in sorted(st, key=str): print(st[k], k)
for k in sorted(pats, key=str)[:40]: print(pats[k], k)
