"""LEAST hypothesis.  Full BFS of states reachable from the census start by single A->B / A->C migrations (cells with |alpha|<=4).
For each edge sigma -> sigma' = sigma;x->T, let G = {tau reachable : nestdata(tau)=nestdata(sigma'), Phi(tau) >= Phi(sigma)}.
Test: sigma' is the unique Phi-minimum of G.  Also report how often |G|>1 (non-vacuity) and whether a same-nestdata pair is Phi-incomparable."""
from qp import *
from collections import Counter, defaultdict
st = Counter(); bad = []
for key, b0 in census():
    if len(b0.nest[0]['pos']) > 4: continue
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    seen = {}; edges = []; stack = [b0]
    while stack:
        b = stack.pop(); k = fullkey(b)
        if k in seen: continue
        seen[k] = b
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); edges.append((k, fullkey(bx), t)); stack.append(bx)
    PH = {k: phi(b) for k, b in seen.items()}
    groups = defaultdict(list)
    for k in seen: groups[k[0]].append(k)
    for g in groups.values():
        st['groups'] += 1; st['groups size>1'] += len(g) > 1
        for i in range(len(g)):
            for j in range(i+1, len(g)):
                st['same-nest pairs'] += 1; st['same-nest pairs comparable'] += leq(PH[g[i]], PH[g[j]]) or leq(PH[g[j]], PH[g[i]])
    for k, k2, t in edges:
        st['edges'] += 1; st[('edge Phi increases', t)] += leq(PH[k], PH[k2])
        G = [m for m in groups[k2[0]] if leq(PH[k], PH[m])]
        st['G size>1'] += len(G) > 1
        mins = [m for m in G if all(leq(PH[m], PH[o]) for o in G)]
        ok = mins == [k2]; st['edge target = unique Phi-min of G'] += ok
        if not ok and len(bad) < 3: bad.append((key, t, len(G)))
    st['cells'] += 1
for k in sorted(st, key=str): print(k, st[k])
print(bad)
