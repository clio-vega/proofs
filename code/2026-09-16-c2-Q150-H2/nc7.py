import sys, random; sys.path.insert(0,'.')
from independent_maya import commutator
from newcounter import tensor_power, Wof
from singleton import crit_star, borders
from itertools import product

def build(d,w,k):
    sig={u:(1 if u==w else 0) for u in product((0,1),repeat=d-1)}
    return Wof(sig,d), Wof(tensor_power(w,d,k),k*d)

def scan_dens(e,W,f,Wb,dens,trials,seed):
    rng=random.Random(seed); width=e+f+4
    for _ in range(trials):
        M=frozenset(i for i in range(width) if rng.random()<dens)
        if commutator({M:1},e,W,f,Wb): return M
    return None

for d in (6,7,8,9):
    bad=[x for x in product((0,1),repeat=d-1) if x[0]!=x[d-2] and not crit_star(x)]
    for w in bad:
        found=None
        for dens in (0.2,0.3,0.4,0.5,0.6,0.7,0.8):
            found=scan_dens(d,*build(d,w,2)[:1],2*d,build(d,w,2)[1],dens,6000,7) if False else None
        # straightforward loop
        W,Wb=build(d,w,2)
        for dens in (0.2,0.3,0.4,0.5,0.6,0.7,0.8):
            found=scan_dens(d,W,2*d,Wb,dens,8000,11)
            if found: break
        print('d=%d w=%s  borders=%s  witness %s' % (d,''.join(map(str,w)),borders(w),
              ('FOUND at density %.1f : M=%s'%(dens,sorted(found))) if found else 'NOT FOUND in 56000 samples'))
        sys.stdout.flush()
