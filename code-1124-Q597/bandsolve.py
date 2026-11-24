"""At the turn vertex vbar (j=e) of split records: the two edges at vbar between proj(tau3) and proj(Z1) (Z1 = tau3 - alpha)
are the TARGET edges.  Certified constraints from bandcert.py (k_i = 0 where certified, k_i >= 1 where a known square
lies over E_i) plus hexagon closure k_i - k_{i+3} = k_{i+4} - k_{i+1} = k_{i+2} - k_{i+5}: are the two target stacks
forced to 0?  (enumeration of k in [0,6]^6).  usage: python3 bandsolve.py SIZE"""
import sys, itertools
from collections import Counter
SIZE = int(sys.argv[1])
src = open('/home/clio/projects/proofs/code-1124-Q597/bandcert.py').read()
src = src.split('for kk in sorted(st')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE)
src = src.replace("        st[(T, 'e' if j == e else 'sh', 'cert', c, 'three-consec-0', three, 'rh', nrh)] += 1",
"""        if j != e: continue
        t3 = next(u for u in Jo[j][1] if zisL(u, k) and v in u[1])
        al = (0, 0, 0, 1) if T == 'B' else (0, -1, 0, 0)
        Z1 = frozenset(tuple(a - b for a, b in zip(Q, al)) for Q in t3[1])
        pt, pz = zpk(t3, k), frozenset(zfree(Q, k) for Q in Z1)
        idx = lambda tri: [i for i, d in enumerate(DIRS) if (vb[0] + d[0], vb[1] + d[1]) in tri]
        it, iz = set(idx(pt)), set(idx(pz))
        tgt = tuple(sorted(it ^ iz)) if len(it & iz) == 0 else None
        # target = the edges strictly between: tau3 spans dirs it, Z1 spans iz; separating triangle shares one dir with each
        sep = sorted(it & iz) if it & iz else None
        between = None
        for i in range(6):
            if {i, (i+1) % 6} == it and {(i+2) % 6, (i+3) % 6} == iz: between = ((i+1) % 6, (i+2) % 6)
            if {i, (i+1) % 6} == iz and {(i+2) % 6, (i+3) % 6} == it: between = ((i+1) % 6, (i+2) % 6)
        forced = None
        if between:
            sols = []
            for K in itertools.product(range(7), repeat=6):
                if any(c[i] == '0' and K[i] for i in range(6)) or any(c[i] == '+' and not K[i] for i in range(6)): continue
                a_, b_, d_ = K[0]-K[3], K[4]-K[1], K[2]-K[5]
                if a_ == b_ == d_: sols.append(K)
            forced = all(K[between[0]] == 0 and K[between[1]] == 0 for K in sols)
        st[(T, 'cert', c, 'tau3 dirs', tuple(sorted(it)), 'Z1 dirs', tuple(sorted(iz)), 'between', between, 'forced0', forced)] += 1""")
exec(src)
for kk in sorted(st, key=str): print(st[kk], kk)
