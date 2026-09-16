import sys
sys.path.insert(0,'.')
from independent_maya import *
from singleton import crit_star
from itertools import product

def tensor_power(w, d, k):
    """sigma = tau^{tensor k}, tau = indicator of w (rank d).  Returns weight dict on {0,1}^{kd-1}."""
    m = k*d
    S = {}
    for y in product((0,1), repeat=m-1):
        val = 1
        for i in range(k):
            blk = y[i*d : i*d + d-1]
            if blk != w: val = 0; break
        S[y] = val
    return S

def Wof(sig, rank):
    return {u: ((-1)**sum(u))*v for u,v in sig.items()}

def check(d, w, k, N):
    e = d; f = k*d
    sig = {u: (1 if u==w else 0) for u in product((0,1),repeat=d-1)}
    sigbar = tensor_power(w, d, k)
    W = Wof(sig,d); Wb = Wof(sigbar,f)
    nz = sum(1 for v in Wb.values() if v)
    bad = 0; acted = 0
    for lam in all_partitions_upto(N):
        S = maya_of_partition(lam)
        if apply_R({S:1}, e, W) or apply_R({S:1}, f, Wb): acted += 1
        if commutator({S:1}, e, W, f, Wb): bad += 1
    return bad, acted, nz

print('=== single-ribbon-shape commuting operators, verified on Maya diagrams ===')
for d,k,N in [(3,2,11),(4,2,12),(5,2,13),(5,3,15),(6,2,14)]:
    good = [w for w in product((0,1),repeat=d-1) if crit_star(w)]
    w = good[0]
    bad, acted, nz = check(d,w,k,N)
    print('d=%d w=%s -> ribbon shape %s ; (e,f)=(%d,%d): |supp Wbar|=%d ; over |lam|<=%d : %d nonzero commutators, %d partitions where some R acts'
          % (d, ''.join(map(str,w)), alpha(w,d), d, k*d, nz, N, bad, acted))
    sys.stdout.flush()

print()
print('=== NEGATIVE CONTROL: singletons that FAIL the border criterion must NOT commute ===')
for d,k,N in [(6,2,14),(7,2,15)]:
    bad_w = [w for w in product((0,1),repeat=d-1) if w[0]!=w[d-2] and not crit_star(w)]
    for w in bad_w[:2]:
        bad, acted, nz = check(d,w,k,N)
        print('d=%d w=%s (w_1!=w_{d-1} but border test fails): %d nonzero commutators  %s'
              % (d, ''.join(map(str,w)), bad, 'FIRES (correct)' if bad else '*** SILENT - PROBLEM ***'))
        sys.stdout.flush()
