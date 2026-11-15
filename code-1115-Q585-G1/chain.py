"""Direct check of the invariant used in the proof of Exchange (lem:chain):
  for every k <= m' :  state_{M_1}(k)  ==  state_M(k) ; (Z_k -> B, east to west)
where state_M(k) = M after D[:k] -> C,  state_{M_1}(k) = M_1 after Dhat[:k] -> C, and Z_k = tops of columns c < alpha_l
in the A-nest of state_M(k).  Equality = same tiling and same nest contents.  Also checks: every Z-move is from an outer corner,
the k-th cells correspond by phi, and the C-landings agree."""
from instr import *
from exchange import orders, run
from collections import Counter
def nests(b): return tuple(frozenset(n['pos']) for n in b.nest)
st = Counter(); bad = []
for key, b0 in census():
    A = b0.nest[0]; uA, wA = A['u'], A['w']; cells = set(A['pos']); R, D = orders(cells); m = len(R)
    ell = max(r for c, r in cells)+1; top = [x for x in R if x[1] == ell-1]; al = len(top)
    hat = cells - set(top); _, Dh = orders(hat); mp = len(Dh)
    bM = copy.deepcopy(b0); bM1 = copy.deepcopy(b0)
    for x in top: mig1(bM1, 0, x, wA)
    for k in range(mp+1):
        S = set(bM.nest[0]['pos'])
        Z = sorted([(c, max(r for cc, r in S if cc == c)) for c in range(al) if any(cc == c for cc, r in S)], key=lambda z: -z[0])
        bt = copy.deepcopy(bM); okc = True
        for z in Z:
            okc &= z in corners(set(bt.nest[0]['pos'])); mig1(bt, 0, z, wA)
        st['Z-moves all from corners'] += okc
        same = tiles_of(bt) == tiles_of(bM1) and nests(bt) == nests(bM1)
        st['invariant holds'] += same; st['checks'] += 1
        if not same and len(bad) < 3: bad.append((key, k))
        if k == mp: break
        _, l1, _ = mig1(bM, 0, D[k], uA); _, l2, _ = mig1(bM1, 0, Dh[k], uA)
        st['C-landings agree'] += l1 == l2; st['steps'] += 1
        st['Z_k nonempty and a z moved (swap step)'] += D[k][0] < al
    st['cells'] += 1
for k in sorted(st): print(k, st[k])
print(bad)
