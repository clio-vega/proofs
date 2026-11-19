"""PROVE 11-19 c2 (xprov2): x-provenance + edge classes relative to x rhombus P(e0,e1); KEEP env = chunks run of read wake tiles, split by case A (-alpha edge of u) / B. usage xprov.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])

def xprov(u, M0, Jo):
    """trace u backwards? We trace forward: tile u of M_P; which step of x's journey Jo (from M0) removes u, and where does it go."""
    if u not in M0: return ('not in M0',)
    for j,(typ,g,n) in enumerate(Jo):
        if u in g:
            vs=[]
            for t in n:
                if t[0]!=u[0]: continue
                v=tuple(bb-aa for aa,bb in zip(cen(u[1]),cen(t[1])))
                if T([u],v)==frozenset([t]): vs.append(v)
            return (typ, j, vs)
    return ('never removed by x',)
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
                    A, B, Jo = (Jy, Jyx, Jx) if side == 0 else (Jx, Jxy, Jy)
                    mA, mB, p, q = middles(A, B)
                    if not bs_test(mA, mB): continue
                    st['records'] += 1
                    MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
                    a, C = mA[0], mA[1:]
                    alpha = disp(a, a); ma = tuple(-v for v in alpha)
                    d = finddir(MP, mA)
                    if d is None: continue
                    # is alpha (up to sign) an edge of x's rhombus at start of Jo?
                    Rx = rhof(Jo[0]); st[(side,'alpha = k*edge of x rhombus', str(mult_of_edge(Rx, alpha)))] += 1
                    st[(side,'y prefix length p', p)] += 1
                    seen=set()
                    for i in range(len(C)):
                        Mi, whoM = made_by(MP, [a] + list(C[:i]))
                        for _, S, _ in candidates(Mi, rhof(C[i]), d):
                            for u in S:
                                if u[0] == 'rh' or u in whoM or u in MxP or u in seen: continue
                                seen.add(u)
                                case = 'A' if mult_of_edge(u, ma) is not None else 'B'
                                xp = xprov(u, M0, Jo)
                                if len(xp)==3:
                                    typ,j,vs = xp
                                    lab = (typ, 'v=-alpha' if ma in vs else ('translated other' if vs else 'not translated'))
                                else: lab = xp
                                st[(side, case, u[0], 'x removes it:', lab)] += 1
                                ecl = classes(u); Rv = [vec(e) for e in tedges(Rx)]
                                epar = [c for c in Rv if par(c, alpha)]; eoth = [c for c in Rv if not par(c, alpha)]
                                st[(side, case, u[0], 'u has edge || alpha', any(par(c, alpha) for c in ecl), '|| other x-edge', any(par(c, eoth[0]) for c in ecl))] += 1
                                if case == 'B' and len(xp) == 3:
                                    st[(side, 'B', u[0], 'x step that moves it: index in run', xp[0], 'x-run vector = -alpha', ma in xp[2])] += 1
                                if case=='A':
                                    # who is u+alpha (behind u)? and u - alpha is common
                                    ub = next(iter(T([u], alpha)))
                                    st[(side,'A', u[0], 'u+alpha in M_P', ub in MP, 'in MxP', ub in MxP)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
