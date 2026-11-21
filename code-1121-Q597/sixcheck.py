"""PROVE 11-21: size-6 fixed subsample (chunks 0..KEEP-1 of NCH, as aw6/pot6) of the (A1) facts:
alpha name, (class,type) of a, Mx step type, first contact == p, contact tiles and Jo's (class,type) at their removal,
y rhombus at a Jo rhombus position in steps 0..p.   usage: KEEP=3 python3 sixcheck.py 6 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
from fractions import Fraction as Fr
UA = (0, -1, 0, 0); WA = (0, 0, 0, 1)
def name(al):
    al = tuple(Fr(x) for x in al)
    for nm, v in (('u_A', UA), ('w_A', WA)):
        if al == tuple(Fr(x) for x in v): return nm
    return 'other'
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
                    st[(tag, 'records')] += 1
                    TA, TB, TO = typed(A), typed(B), typed(Jo)
                    st[(tag, 'alpha', name(disp(A[p], A[p])), 'a', TA[p], 'Mx', TB[p])] += 1
                    W = set(M0) - set(apply(M0, Jo))
                    fc = next((k for k, s in enumerate(A) if any(u in W for u in s[1] if u[0] != 'rh')), None)
                    st[(tag, 'first contact == p', fc == p)] += 1
                    rem = {}
                    for j, s in enumerate(Jo):
                        for u in s[1]: rem.setdefault(u, j)
                    st[(tag, 'contact tiles with Jo config', tuple(sorted((u[0],) + TO[rem[u]] for u in A[p][1] if u in W)))] += 1
                    JR = {frozenset(rh(s[1])) for s in Jo} | {frozenset(rh(Jo[-1][2]))}
                    st[(tag, 'y rhombus at a Jo rhombus position, steps 0..p', any(frozenset(rh(s[1])) in JR for s in A[:p+1]))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns, 'chunks', os.environ.get('KEEP', NCH), 'of', NCH)
    for k in sorted(st, key=str): print(st[k], k)
