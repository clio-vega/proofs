"""Value picture on I: the set of vertex values carried by I-cells in M, Mx, My, Mxy, relative to min value in M."""
import pickle
from collections import Counter
exec(open('chains.py').read().split("BC = [r for r")[0])
st = Counter()
for r in recs:
    if r['kind'] != 'BC': continue
    X, Y, I = solo(r)
    if not I: continue
    side = next(iter(I))[0]
    Xv, Yv = vec(r['X']), vec(r['Y'])
    val = lambda k, *vs: tuple(r['P'][k][i] + sum(v.get(k, Z)[i] for v in vs) for i in (0, 1))
    s = {0: (0, -1), 1: (0, 1)}[side]
    VM = {r['P'][k] for k in I}
    VMxy = {tuple(r['P'][k][i] + s[i] for i in (0, 1)) for k in I}
    VMx = {val(k, Xv) for k in I}; VMy = {val(k, Yv) for k in I}
    base = min(VM)
    rel = lambda S: tuple(sorted((a-base[0], b-base[1]) for a, b in S))
    st[('q' if side == 0 else 'p', 'M', rel(VM), 'Mx', rel(VMx), 'My', rel(VMy), 'Mxy', rel(VMxy))] += 1
for k, v in sorted(st.items(), key=lambda kv: -kv[1]): print(v, k)
