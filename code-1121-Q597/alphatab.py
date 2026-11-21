"""PROVE 11-21 Step 0 (Q597 G1a, piece A1): what is alpha?
alpha = rhombus displacement of the first MIDDLE step a of the clean journey (same recipe as g1aparity.py / bs6.py).
Per record: side, clean target, prefix length p, prefix word, type of a, class of R_a, alpha, alpha in {u_A, w_A},
-alpha edge vector of the OTHER journey's R_0, nest frame (u,w) as the engine stores it, number of A-tiles.
Preconditions asserted: record clean nontrivial (bs_test), len(middle) >= 2.
Control: B/B non-commuting pairs, both journeys, same recipe.
usage: python3 alphatab.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
from fractions import Fraction as Fr
UA = (0, -1, 0, 0); WA = (0, 0, 0, 1)
def nv(v): return tuple(Fr(x) for x in v)
def edgevec(t, v): return mult_of_edge(t, v) is not None
def rclass(s):
    """edge vectors of the rhombus before step s (the 'g' rhombus), as a sorted tuple of Lambda_0 part"""
    vs = sorted(rh(s[1])); P = vs[0]
    ds = sorted({tuple(Fr(b) - Fr(a) for a, b in zip(P, Q)) for Q in vs if Q != P})
    return ds
def fr(v): return tuple(str(x) for x in v)
def name(al):
    al = nv(al)
    if al == nv(UA): return 'u_A'
    if al == nv(WA): return 'w_A'
    if al == tuple(-x for x in nv(UA)): return '-u_A'
    if al == tuple(-x for x in nv(WA)): return '-w_A'
    return 'other' + str(fr(al))
def word(J): return ''.join('b' if s[0] == 'b' else ('a' + str(s[0])[-1] if s[0] != 'a' else 'a') for s in J)
def work(w):
    st = Counter(); nstates = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            nstates += 1
            st[('frame u,w', fr(b.nest[0]['u']), fr(b.nest[0]['w']))] += 1
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        if kind == 'BB':
                            if tiles_of(b12) == tiles_of(b21): continue
                            for nm, A_, B_, Jo_ in (('y', Jy, Jyx, Jx), ('x', Jx, Jxy, Jy)):
                                mA_, mB_, p_, q_ = middles(A_, B_)
                                if len(mA_) < 2: st[('BB', 'middle<2')] += 1; continue
                                al = disp(mA_[0], mA_[0])
                                st[('BB', nm, 'p', p_, 'type a', mA_[0][0], 'alpha', name(al), '-alpha edge of Jo R0', edgevec(rhof(Jo_[0]), tuple(-v for v in al)))] += 1
                            continue
                        P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                        X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                        K = set(X) | set(Y) | set(Yx) | set(Xy)
                        I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                        sides = {k[0] for k in I}
                        if len(sides) != 1: continue
                        side = sides.pop()
                        A, B, Jo = (Jy, Jyx, Jx) if side == 0 else (Jx, Jxy, Jy)
                        Tc = 'C' if side == 0 else 'B'
                        mA, mB, p, q = middles(A, B)
                        if not bs_test(mA, mB): continue
                        assert len(mA) >= 2
                        st['records'] += 1
                        a = mA[0]; al = disp(a, a)
                        ra = rclass(a); r0c = rclass(A[0]); r0o = rclass(Jo[0])
                        st[(side, 'clean->' + Tc, 'p', p, 'prefix', word(A[:p]), 'type a', a[0], 'alpha', name(al))] += 1
                        st[(side, 'alpha', name(al), 'R_a same shape as clean R_0', ra == r0c, 'clean R_0 same shape as Jo R_0', r0c == r0o)] += 1
                        st[(side, 'alpha', name(al), '-alpha edge of Jo R_0', edgevec(rhof(Jo[0]), tuple(-v for v in al)), '-alpha edge of clean R_0', edgevec(rhof(A[0]), tuple(-v for v in al)), 'edge of R_a', edgevec(rhof(a), tuple(-v for v in al)))] += 1
                        st[(side, 'first type of clean', A[0][0], 'first type of Jo', Jo[0][0], 'alpha', name(al))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
