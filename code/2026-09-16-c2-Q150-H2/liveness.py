import sys; sys.path.insert(0,'.')
from lib import *
from independent_maya import apply_R
print('=== exhaustive window scan: commutator zero AND two-bead sector genuinely live ===')
for d,k,width in [(3,2,13),(4,2,15)]:
    for w in [x for x in product((0,1),repeat=d-1) if crit_star(x)]:
        W,Wb=build(d,w,k); e=d; f=k*d
        bad=0; live=0; tot=0
        for bits in product((0,1),repeat=width):
            M=frozenset(i for i in range(width) if bits[i])
            tot+=1
            a=apply_R(apply_R({M:1},e,W),f,Wb)
            b=apply_R(apply_R({M:1},f,Wb),e,W)
            # liveness: does EITHER order reach a state differing from M in >=4 sites?
            if any(len(set(T)^set(M))>=4 for T in list(a)+list(b)): live+=1
            out={}
            for kk,v in a.items(): out[kk]=out.get(kk,0)+v
            for kk,v in b.items(): out[kk]=out.get(kk,0)-v
            if any(v for v in out.values()): bad+=1
        print('  d=%d k=%d w=%s : %d/%d nonzero commutators ; %d states reach the two-bead sector'
              % (d,k,''.join(map(str,w)),bad,tot,live))
        sys.stdout.flush()
