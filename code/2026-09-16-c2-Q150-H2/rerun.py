import sys; sys.path.insert(0,'.')
from mayatest import scan
from newcounter import tensor_power, Wof
from singleton import crit_star
from itertools import product

def build(d, w, k):
    sig = {u:(1 if u==w else 0) for u in product((0,1),repeat=d-1)}
    return Wof(sig,d), Wof(tensor_power(w,d,k), k*d)

print('=== POSITIVE: border-criterion singletons  (must be 0) ===')
for d,k in [(3,2),(4,2),(5,2),(6,2),(7,2),(5,3)]:
    for w in [x for x in product((0,1),repeat=d-1) if crit_star(x)][:3]:
        W,Wb = build(d,w,k)
        nz,tb = scan(d,W,k*d,Wb,trials=3000,seed=1)
        print('  d=%d k=%d w=%s : %d/3000 Maya states with nonzero commutator' % (d,k,''.join(map(str,w)),nz))
    sys.stdout.flush()
print()
print('=== NEGATIVE CONTROL: singletons failing the border criterion (must FIRE) ===')
for d,k in [(6,2),(7,2)]:
    for w in [x for x in product((0,1),repeat=d-1) if x[0]!=x[d-2] and not crit_star(x)]:
        W,Wb = build(d,w,k)
        nz,tb = scan(d,W,k*d,Wb,trials=3000,seed=1)
        print('  d=%d w=%s : %d/3000 nonzero (%d involve a two-bead move)  %s'
              % (d,''.join(map(str,w)),nz,tb,'FIRES' if nz else '*** SILENT ***'))
    sys.stdout.flush()
print()
print('=== NEGATIVE CONTROL 2: w with w_1 == w_{d-1} (must FIRE) ===')
for d,k in [(4,2),(5,2)]:
    for w in [x for x in product((0,1),repeat=d-1) if x[0]==x[d-2]][:3]:
        W,Wb = build(d,w,k)
        nz,tb = scan(d,W,k*d,Wb,trials=3000,seed=1)
        print('  d=%d w=%s : %d/3000 nonzero  %s' % (d,''.join(map(str,w)),nz,'FIRES' if nz else '*** SILENT ***'))
    sys.stdout.flush()
