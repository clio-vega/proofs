"""Does (p,q) jointly determine the fine mixed subdivision, over ALL subdivisions?  Uses all_a_b.pkl (mixsub enumeration)."""
import sys, pickle
from collections import defaultdict, Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
for fn in sys.argv[1:]:
    a, b, M = pickle.load(open(fn, 'rb'))
    g = defaultdict(set)
    for T in M:
        key = frozenset(vs for k, vs in T if k in ('t0', 't1'))
        g[key].add(T)
    print(fn, 'subdivisions', len(M), 'pure-cell classes', len(g), 'collisions', sum(len(s) > 1 for s in g.values()),
          Counter(len(s) for s in g.values()))
