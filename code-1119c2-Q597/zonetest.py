"""PROVE 11-19 c2 Step 0: the de Bruijn zone falsifier for WS-b.   usage: python3 zonetest.py {5|6} NW NCH
Zone (de Bruijn), defined on a bare board (set of (kind, frozenset(vertices))):
  a tile with a pair of PARALLEL opposite edges (sq, rh) is crossed from one edge of class c to the other edge of class c;
  the zone Z_c(t) is the maximal chain so obtained; it stops at a tile with no second c-edge (every tri) or at the boundary.
  Triangles lie in NO zone (their three edges are pairwise non-parallel).
For every clean nontrivial B/C record and every read WAKE tile u (non-rh tile of a candidate of y's C-step in M_i, u in M_P \\ Mx_P):
  kind of u; is alpha (or a multiple) an edge vector of u; the translate chain u - j alpha (j=0..) membership in M_P / Mx_P;
  adjacency of u and u - alpha; for squares: zone membership of the chain in M_P and in Mx_P.
Control: B/B journeys, the first read tile whose -alpha translate is ABSENT (the WS-b failure).
Lesion: every square of M_P with alpha replaced by each other edge vector of the square / by alpha itself on non-read squares."""
import sys
exec(open('/home/clio/projects/proofs/code-1119-Q597/prov.py').read().split('\ndef work(w):')[0])
from stepset import tedges, index
def vec(e):
    P, Q = tuple(e); return tuple(Fr(b) - Fr(a) for a, b in zip(P, Q))
def par(v, w):
    """v = k w for some nonzero rational k"""
    k = migrate.Board.coef(v, w) if any(w) else None
    return k is not None and k != 0
def classes(t):
    out = []
    for e in tedges(t):
        v = vec(e)
        if not any(par(v, w) for w in out): out.append(v)
    return out
def cedges(t, c): return [e for e in tedges(t) if par(vec(e), c)]
def zone(Bd, E, t, c):
    """ordered list of tiles of the de Bruijn zone through t crossing edges of class c (t must have two c-edges)"""
    es = cedges(t, c)
    if len(es) != 2: return None
    Z = [t]
    for e0 in es:
        cur, e = t, e0; side = []
        for _ in range(200):
            nb = [s for s in E.get(e, ()) if s != cur]
            if not nb: break
            nxt = nb[0]; ce = cedges(nxt, c)
            side.append(nxt)
            if len(ce) != 2: break             # triangle (or degenerate): zone ends here
            e = ce[0] if ce[1] == e else ce[1]; cur = nxt
        Z = (list(reversed(side)) + Z) if e0 == es[0] else (Z + side)
    return Z
def endtile(Z): return Z[0][0] != 'tri' and Z[-1][0] != 'tri'
def mult_of_edge(t, v):
    for e in tedges(t):
        w = vec(e); k = migrate.Board.coef(v, w)
        if k is not None and k != 0: return k
    return None
def sharevert(s, t): return len(set(s[1]) & set(t[1]))
def chain_sig(u, mal, MP, MxP, jmax=8):
    sig = []
    for j in range(jmax):
        t = next(iter(T([u], tuple(j*x for x in mal))))
        a, b = t in MP, t in MxP
        if not a and not b: break
        sig.append('C' if a and b else ('W' if a else 'X'))
    return ''.join(sig)
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
                    if x < y:   # B/B control (same recipe as prov.py)
                        c1, Kx = mvl(b, x, 'B', DV); c2, Ky = mvl(b, y, 'B', DV)
                        c12, Kyx = mvl(c1, y, 'B', DV); c21, Kxy = mvl(c2, x, 'B', DV)
                        if tiles_of(c12) != tiles_of(c21):
                            for nm, A_, B_, Jo_ in (('y', Ky, Kyx, Kx), ('x', Kx, Kxy, Ky)):
                                mA_, mB_, p_, q_ = middles(A_, B_)
                                if len(mA_) < 2: continue
                                MP_ = apply(M0, A_[:p_]); MxP_ = apply(apply(M0, Jo_), B_[:p_])
                                d_ = finddir(MP_, mA_)
                                if d_ is None: continue
                                al = disp(mA_[0], mA_[0]); mal = tuple(-v for v in al)
                                done = False
                                for i in range(min(len(mA_)-1, len(mB_))):
                                    Mi, whoM = made_by(MP_, list(mA_[:i+1])); Xi, _ = made_by(MxP_, list(mB_[:i]))
                                    for _, S, _ in candidates(Mi, rhof(mA_[i+1]), d_):
                                        for u in S:
                                            if u[0] == 'rh' or u in whoM: continue
                                            u2 = next(iter(T([u], mal)))
                                            if u2 in Xi: continue
                                            # the WS-b failure tile
                                            k = mult_of_edge(u, mal)
                                            st[('BB', nm, 'absent-translate tile', u[0], 'alpha mult of edge', str(k), 'chain', chain_sig(u, mal, MP_, MxP_))] += 1
                                            done = True; break
                                        if done: break
                                    if done: break
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
                    EP, EX = index(MP), index(MxP)
                    a, C = mA[0], mA[1:]
                    alpha = disp(a, a); ma = tuple(-v for v in alpha)
                    d = finddir(MP, mA)
                    if d is None: continue
                    st[('alpha is a multiple of an edge of the a-rhombus', str(mult_of_edge(rhof(a), alpha)))] += 1
                    seen = set()
                    for i in range(len(C)):
                        Mi, whoM = made_by(MP, [a] + list(C[:i]))
                        for _, S, _ in candidates(Mi, rhof(C[i]), d):
                            for u in S:
                                if u[0] == 'rh' or u in whoM or u in MxP: continue    # read WAKE tiles only
                                if u in seen: continue
                                seen.add(u)
                                k = mult_of_edge(u, ma)
                                sig = chain_sig(u, ma, MP, MxP)
                                u1 = next(iter(T([u], ma)))
                                st[(side, 'wake', u[0], '-alpha = k*edge, k', str(k))] += 1
                                st[(side, 'wake', u[0], 'chain', sig)] += 1
                                st[(side, 'wake', u[0], 'shared vertices with u-alpha', sharevert(u, u1))] += 1
                                if u[0] == 'sq':
                                    # zone of M_P through u in the class transverse to alpha (crossing edges NOT parallel to alpha)
                                    cl = classes(u)
                                    tz = [c for c in cl if not par(ma, c)]
                                    for c in tz:
                                        ZP = zone(MP, EP, u, c)
                                        inP = [j for j in range(len(sig)) if sig[j] in 'WC']
                                        tiles_j = [next(iter(T([u], tuple(j*x for x in ma)))) for j in range(len(sig))]
                                        okP = all(tiles_j[j] in ZP for j in inP)
                                        okX = None
                                        if sig[1:2] in ('C', 'X'):
                                            ZX = zone(MxP, EX, tiles_j[1], c)
                                            inX = [j for j in range(len(sig)) if sig[j] in 'XC']
                                            okX = ZX is not None and all(tiles_j[j] in ZX for j in inX)
                                        st[(side, 'wake sq', 'transverse class: chain in one zone of M_P', okP, 'of Mx_P', okX)] += 1
                                        st[(side, 'wake sq', 'zone of M_P length', len(ZP), 'ends', ZP[0][0], ZP[-1][0])] += 1
                    # lesion: every square of M_P not read, translate chain by -alpha: how often is it a zone strip?
                    for u in MP:
                        if u[0] != 'sq' or u in seen: continue
                        k = mult_of_edge(u, ma)
                        u1 = next(iter(T([u], ma)))
                        st[('lesion: unread sq of M_P', '-alpha = k*edge', str(k), 'u-alpha in M_P', u1 in MP)] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
