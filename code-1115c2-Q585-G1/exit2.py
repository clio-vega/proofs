"""Split the clean-exit statistics by WHICH move is second: B-move crossing a C-wake (role 'B after C') or C-move crossing a
B-wake ('C after B').  'rejoin' = nonempty common suffix with the solo journey (as hexagon tile sets)."""
import pickle
from collections import Counter
D = pickle.load(open('/tmp/c2-diamonds.pkl','rb')); st = Counter()
for d in D:
    for second, alone, role in (('h1y','h2y', d['kind'][1]+' after '+d['kind'][0]), ('h2x','h1x', d['kind'][0]+' after '+d['kind'][1])):
        H2 = [s[1] for s in d[second]]; HA = [s[1] for s in d[alone]]
        if H2 == HA: continue
        k = 0
        while k < min(len(H2), len(HA)) and H2[-1-k] == HA[-1-k]: k += 1
        st[(d['kind'], role, 'changed')] += 1; st[(d['kind'], role, 'nonempty common suffix')] += k > 0
for x in sorted(st): print(x, st[x])
