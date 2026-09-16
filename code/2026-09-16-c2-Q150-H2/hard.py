import sys; sys.path.insert(0,'.')
from lib import *
print('=== HARD negative controls: w_1 != w_{d-1} yet the border criterion fails ===')
tot=0; fired=0
for d in range(4,11):
    hard=[x for x in product((0,1),repeat=d-1) if x[0]!=x[d-2] and not crit_star(x)]
    for w in hard:
        ds=[dl for dl in fail_deltas(w) if 1<=dl<=d-1]
        delta=ds[0]; W,Wb=build(d,w,2)
        M,y=witness_state(w,delta,2)
        c=commutator({M:1},d,W,2*d,Wb)
        tb=any(len(set(T)^set(M))>=4 for T in c)
        tot+=1; fired += 1 if c else 0
        print('  d=%2d w=%-10s borders=%-10s delta=%d : %-22s two-bead=%s'
              % (d,''.join(map(str,w)),borders(w),delta,
                 'NONZERO '+str(sorted(c.values())) if c else '*** ZERO - PROBLEM ***', tb))
    sys.stdout.flush()
print('  -> %d/%d hard controls fired' % (fired,tot))
