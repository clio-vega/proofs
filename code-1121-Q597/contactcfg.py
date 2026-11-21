"""PROVE 11-21: contact configurations. At the first contact (y's step a at R_a, which removes u in W), record for each
contact tile u: kind, y's (class,type), Jo's (class,type) at the step j_u removing u, and the offset base(R_a) - base(R_{j_u})
(bases = coordinatewise-min vertex), plus y's previous step (class,type).   usage: python3 contactcfg.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def base(s):
    vs = rh(s[1]); return tuple(min(Fr(v[i]) for v in vs) for i in range(4))
def work(w):
    st = Counter(); nstates = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            nstates += 1; M0 = bare(b)
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    b1, Jx = mvl(b, x, 'B', DV); b2, Jy = mvl(b, y, 'C', DV)
                    b12, Jyx = mvl(b1, y, 'C', DV); b21, Jxy = mvl(b2, x, 'B', DV)
                    P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                    X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                    K = set(X) | set(Y) | set(Yx) | set(Xy)
                    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                    sides = {k[0] for k in I}
                    if len(sides) != 1: continue
                    side = sides.pop()
                    tag, A, B, Jo = (('clean->C', Jy, Jyx, Jx) if side == 0 else ('clean->B', Jx, Jxy, Jy))
                    mA, mB, p, q = middles(A, B)
                    if not bs_test(mA, mB): continue
                    W = set(M0) - set(apply(M0, Jo))
                    TA = typed(A); TO = typed(Jo)
                    rem = {}
                    for j, s in enumerate(Jo):
                        for u in s[1]: rem.setdefault(u, j)
                    for u in A[p][1]:
                        if u not in W: continue
                        j = rem[u]
                        off = tuple(int(c - d) for c, d in zip(base(A[p]), base(Jo[j])))
                        st[(tag, u[0], 'y', TA[p], 'Jo', TO[j], 'off', off)] += 1
                    st[(tag, 'prev y step', TA[p-1] if p else None)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
