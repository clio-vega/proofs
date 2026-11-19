"""Is -alpha the zone vector of x's runs that touch y's crossing?  For each clean record: x's type-a steps whose tiles meet
the tiles consumed by y's M-middle; collect the displacement vectors of their displaced non-rh tiles, split by zone."""
exec(open('runzone.py').read().split('st = Counter(); seen = set()')[0])
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, A, B, Jo = clean(r); mA, mB, p, q = middles(A, B)
    a = mA[0]; alpha = tuple(y-x for x, y in zip(cen(rh(a[1])), cen(rh(a[2])))); ma = tuple(-v for v in alpha)
    used = set().union(*[s[1] for s in mA]) | set().union(*[s[2] for s in mA])
    zv = Counter()
    for s in Jo:
        if s[0] != 'a': continue
        e = ev(('rh', rh(s[1]))); e0 = next(v for v in e if v[2:] == (0, 0)); e1 = next(v for v in e if v[:2] == (0, 0))
        for t in s[1]:
            if t[0] == 'rh' or t not in used: continue
            vs = [v for v in {tuple(b-a_ for a_, b in zip(min(t[1]), min(u[1]))) for u in s[2] if u[0] == t[0]} if Tt(t, v) in s[2]]
            if vs: zv[('e1zone' if e1 in ev(t) else 'e0zone', vs[0] == ma)] += 1
    st[(side, tuple(sorted(zv)))] += 1
for k in sorted(st, key=str): print(st[k], k)
