"""Under BS: contribution of each crossing step s in C is T_{-alpha} delta_s - delta_s; of a: T_gamma delta_a - delta_a.
Check: (i) alpha is supported on one coordinate block (P0 on q-side, P1 on p-side); (ii) steps whose delta lives on the OTHER side
contribute 0 (automatic: T shifts their keys by the zero block); (iii) delta_a lives on the interaction side;
(iv) number of crossing steps with nonzero contribution."""
from collections import Counter
from q597lib import *
exec(open('reduce.py').read().split('st = Counter()')[0].split('from q597lib import *')[1])
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side = r['Us'][0][0]
    A, B = (r['Jy'], r['Jyx']) if side == 0 else (r['Jx'], r['Jxy'])
    mA, mB, p, q = middles(A, B); a, C = mA[0], mA[1:]
    alpha = tuple(b-x for x, b in zip(cen(rh(a[1])), cen(rh(a[2]))))
    gamma = tuple(b-x for x, b in zip(cen(rh(C[0][1])), cen(rh(C[-1][2]))))
    blk = 'P0' if any(alpha[:2]) and not any(alpha[2:]) else ('P1' if any(alpha[2:]) and not any(alpha[:2]) else 'mixed')
    st[(side, 'alpha block', blk)] += 1
    st[(side, 'delta_a sides', tuple(sorted({k[0] for k in sdelta(a)})))] += 1
    nz = 0
    for s in C:
        c = addd(T(sdelta(s), tuple(-x for x in alpha)), neg(sdelta(s)))
        sides = tuple(sorted({k[0] for k in sdelta(s)}))
        st[(side, 'C-step delta sides', sides, 'contrib zero', not c)] += 1
        nz += bool(c)
    st[(side, '#C steps nonzero', nz, 'of', len(C))] += 1
for k in sorted(st, key=str): print(st[k], k)
