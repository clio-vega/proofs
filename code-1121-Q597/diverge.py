"""PROVE 11-21 Step 0b: the divergence step a of the clean journey (first middle step in M) versus the step c_1 - alpha
taken from the same rhombus R_a in Mx.  Typed words (class index 0 = start class, 1 = other; type a0/a1/b).
Per record: typed prefix, (class,type) of a, (class,type) of Mx's step p, is a the FIRST (class,type) of its kind in the journey,
which tiles of a's removed set g_a are absent from Mx_P, and is c_1-alpha's configuration available in M_P (both-available => tie-break).
Control: B/B both journeys.   usage: python3 diverge.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
from fractions import Fraction as Fr
def l0(s):
    """Lambda_0 edge direction of the rhombus removed by step s (the rhombus before s), up to sign"""
    vs = list(rh(s[1]))
    for P in vs:
        for Q in vs:
            d = tuple(Fr(b) - Fr(a) for a, b in zip(P, Q))
            if any(d[:2]) and not any(d[2:]):
                return d[:2] if (d[0], d[1]) > (0, 0) else tuple(-x for x in d[:2])
def typed(J):
    c0 = l0(J[0]); out = []
    for s in J:
        cl = 0 if l0(s) == c0 else 1
        if s[0] == 'b': t = 'b'
        else:
            al = disp(s, s); t = 'a1' if not any(al[:2]) else ('a0' if not any(al[2:]) else 'a?')
        out.append((cl, t))
    return out
def ws(T): return ' '.join('%d%s' % ct for ct in T)
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
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        if kind == 'BB':
                            if tiles_of(b12) == tiles_of(b21): continue
                            roles = (('BBy', Jy, Jyx, Jx), ('BBx', Jx, Jxy, Jy))
                        else:
                            P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                            X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                            K = set(X) | set(Y) | set(Yx) | set(Xy)
                            I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                            sides = {k[0] for k in I}
                            if len(sides) != 1: continue
                            side = sides.pop()
                            roles = (((side, 'clean->C'), Jy, Jyx, Jx),) if side == 0 else (((side, 'clean->B'), Jx, Jxy, Jy),)
                        for tag, A, B, Jo in roles:
                            mA, mB, p, q = middles(A, B)
                            if kind == 'BC' and not bs_test(mA, mB): continue
                            if len(mA) < 2 or len(mB) < 1: st[(tag, 'short middle')] += 1; continue
                            if kind == 'BC': st[(tag, 'records')] += 1
                            TA = typed(A); TB = typed(B)
                            ta, tb = TA[p], TB[p]
                            first = TA.index(ta) == p
                            st[(tag, 'a', ta, 'Mx step p', tb, 'a is first of its (class,type)', first)] += 1
                            st[(tag, 'prefix', ws(TA[:p]), 'a', ta, 'Mx', tb)] += 1
                            MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
                            assert rh(A[p][1]) == rh(B[p][1])
                            miss = [u for u in A[p][1] if u[0] != 'rh' and u not in MxP]
                            gb_in_M = all(u in MP for u in B[p][1] if u[0] != 'rh')
                            st[(tag, 'g_a tiles absent from Mx_P', len(miss), 'Mx config available in M_P', gb_in_M)] += 1
                            # where were the missing tiles removed by Jo: step index, class parity, type
                            for u in miss:
                                jr = next((j for j, s in enumerate(Jo) if u in s[1]), None)
                                TO = typed(Jo)
                                st[(tag, 'missing tile', u[0], 'removed by Jo at', None if jr is None else (TO[jr], 'b-parity', sum(1 for s in Jo[:jr] if s[0] == 'b') % 2, 'j==0', jr == 0))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
