"""Window rewrite rule: for each changed journey, the solo window word vs the crossed window word (types a/b, after stripping the
longest common prefix and suffix of the HEXAGON sequence).  Tests: reverse, same multiset, cyclic rotation."""
import pickle
from collections import Counter
D = pickle.load(open('/tmp/c2-diamonds.pkl','rb')); st = Counter(); pairs = Counter()
for d in D:
    if d['kind'] != 'BC': continue
    for second, alone, role in (('h1y','h2y','C after B'), ('h2x','h1x','B after C')):
        H2, HA = d[second], d[alone]
        if [s[1] for s in H2] == [s[1] for s in HA]: continue
        p = 0
        while p < min(len(H2), len(HA)) and H2[p][1] == HA[p][1]: p += 1
        q = 0
        while q < min(len(H2), len(HA)) - p and H2[-1-q][1] == HA[-1-q][1]: q += 1
        w2 = ''.join(s[0] for s in H2[p:len(H2)-q]); wa = ''.join(s[0] for s in HA[p:len(HA)-q])
        st[(role, 'n')] += 1
        st[(role, 'crossed = reverse(solo)')] += w2 == wa[::-1]
        st[(role, 'same multiset')] += sorted(w2) == sorted(wa)
        st[(role, 'cyclic rotation')] += len(w2) == len(wa) and w2 in wa + wa
        pairs[(role, wa, w2)] += 1
for k in sorted(st): print(k, st[k])
for k, v in pairs.most_common(14): print(v, k)
