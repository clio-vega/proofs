"""Triangle slide: in tiling T, a 0-triangle tau x {w} and the parallelograms e x f (e an edge of tau, f=[w,w']) on one side of
tau (w.r.t. f) together tile the pentagon tau+f; replace by tau x {w'} and the parallelograms on the other side.  p is unchanged,
q(tau) changes w -> w'.  Every candidate is re-validated by control_pool.valid (independent of how it was produced)."""
import sys, pickle, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E, control_pool as CP
from collections import Counter
E.SIGN = 1; f0 = f1 = CP.MAPS['(i+j,i)']
def slides(cells, b):
    """cells: set of frozensets (my coordinates).  yield flipped cell sets."""
    E1 = E.unit_edges(b)
    for c in cells:
        B0, B1 = E.parts(c)
        if len(B1) != 1 or len(B0) != 3: continue
        w = next(iter(B1))
        for f in E1:
            if w not in f: continue
            (w2,) = tuple(f - {w})
            edges = [frozenset(e) for e in itertools.combinations(sorted(B0), 2)]
            par = {e: frozenset(x + y for x in e for y in f) for e in edges}
            present = [e for e in edges if par[e] in cells]
            if not present: continue
            new_tri = frozenset(x + w2 for x in B0)
            rest = [e for e in edges if e not in present]
            out = set(cells) - {c} - {par[e] for e in present} | {new_tri} | {par[e] for e in rest}
            yield frozenset(out), (B0, w, w2, len(present))
pool = pickle.load(open('pool.pkl', 'rb'))
st = Counter(); example = None
for (a, b), D in sorted(pool.items()):
    for T in D:
        cells = set(frozenset(CP.conv(v, f0, f1) for v in vs) for _, vs in T)
        kinds = {vs: k for k, vs in E.make_tiles(a, b)}
        found = 0
        for new, info in slides(cells, b):
            tl = [(kinds.get(c, '?'), tuple(sorted(c))) for c in new]
            if any(k == '?' for k, _ in tl): continue
            # validate in my coordinates: build a pseudo-engine tiling and use valid with identity maps
            ok = CP.valid([(k, vs) for k, vs in tl], a, b, lambda m, n: (m, n), lambda m, n: (m, n))
            st[((a, b), ok)] += 1
            if ok == 'ok':
                found += 1
                if example is None: example = ((a, b), sorted(map(sorted, cells)), sorted(map(sorted, new)), info)
        st[((a, b), 'tilings with >=1 valid slide' if found else 'tilings with no valid slide')] += 1
    if a + b >= 10: break
for k in sorted(st, key=str): print(k, st[k])
pickle.dump(example, open('slide_example.pkl', 'wb'))
print('example info', example[0], example[3] if example else None)
