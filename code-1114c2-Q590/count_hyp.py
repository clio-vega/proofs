"""Hypothesis K: over ALL full mixed emptyings of A (each step: any current corner, to B or to C), the B-landing set L_B
depends only on |S_B| (how many cells went to B).  Weaker K': depends only on the SET S_B.  K for all k => Order Reversal
(compare 'R[:m-k] to B, rest to C' with 'D[:k] to C, rest to B').  Also test the specific pair directly (OR-k)."""
from common import *
from collections import Counter, defaultdict
st = Counter(); wit = []
def key_of(b): return (frozenset(b.nest[0]['pos']), tuple(tiles_of(b)))
for key, b0 in census():
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}; beta0 = set(b0.nest[1]['pos'])
    # DFS over states carrying (S_B frozenset, L_B frozenset); memo on (state, S_B, L_B)
    res_k = defaultdict(set); res_S = defaultdict(set); seen = set()
    stack = [(b0, frozenset())]
    while stack:
        b, SB = stack.pop(); k = (key_of(b), SB)
        if k in seen: continue
        seen.add(k); cs = corners(set(b.nest[0]['pos']))
        if not cs:
            LB = frozenset(set(b.nest[1]['pos']) - beta0); res_k[len(SB)].add(LB); res_S[SB].add(LB); continue
        for x in cs:
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append((bx, SB | {x} if t == 'B' else SB))
    st['cells'] += 1
    kfail = sum(len(v) > 1 for v in res_k.values()); Sfail = sum(len(v) > 1 for v in res_S.values())
    st['K: (cell,k) classes'] += len(res_k); st['K: classes with >1 L_B'] += kfail
    st["K': (cell,S_B) classes"] += len(res_S); st["K': classes with >1 L_B"] += Sfail
    st['cells K holds'] += kfail == 0; st["cells K' holds"] += Sfail == 0
    if kfail and len(wit) < 2: wit.append((key, {k: len(v) for k, v in res_k.items()}))
print(dict(st)); print(wit)
