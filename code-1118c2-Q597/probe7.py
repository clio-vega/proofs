"""Shared track: is the run w (C minus its two b-turns) of y's M-middle a rigid translate (one vector delta) of a contiguous block of x's journey Jo?
Also: same question for y's Mx-middle run vs x's journey in M (vector delta - alpha?)."""
from replay import *
from collections import Counter
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, A, B, Jo = clean(r); mA, mB, p, q = middles(A, B)
    w = mA[2:-1]; wx = mB[1:-2]
    if not w: st[(side, 'empty run')] += 1; continue
    hits = []
    for i in range(len(Jo)-len(w)+1):
        v = tvec(Jo[i], w[0])
        if v is not None and all(tr(Jo[i+j], v) == w[j] for j in range(len(w))): hits.append((i, v))
    st[(side, 'run = translate of x-block', bool(hits))] += 1
    if hits:
        st[(side, 'delta', tuple(str(c) for c in hits[0][1]))] += 1
        # does the x-block extend (one more step either side) as translate?
        i, v = hits[0]; L = len(w)
        st[(side, 'x-block extends before', i > 0 and tr(Jo[i-1], v) == mA[1])] += 1
for k in sorted(st, key=str): print(st[k], k)
