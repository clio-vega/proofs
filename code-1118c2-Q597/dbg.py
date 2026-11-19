from replay import *
from stepset import *
import cmath
r = next(r for r in R if r['kind'] == 'BC' and r['Us'])
side, A, B, Jo = clean(r); mA, mB, p, q = middles(A, B)
MP = apply(board0(r), A[:p]); s = mA[0]; Rm = next(u for u in s[1] if u[0] == 'rh')
print('recorded', s[0], len(s[1]), Rm in MP)
for k in range(12):
    d = cmath.exp(1j*cmath.pi*k/6); c = candidates(MP, Rm, d)
    got, _ = step(MP, Rm, d)
    print(k, len(c), [x[0] for x in c], got is not None and got[1] == s[1], got is not None and got[2] == s[2])
print('---')
d = cmath.exp(0j); Bm = set(MP)
for i, s in enumerate(mA):
    Rm = next(u for u in s[1] if u[0] == 'rh'); got, c = step(Bm, Rm, d)
    print(i, s[0], [x[0] for x in c], None if got is None else (got[1] == s[1], got[2] == s[2]))
    Bm = apply(Bm, [s])
print([ (k, vs) for k, vs in r['M'].values() if k == 'rh'][:2])
