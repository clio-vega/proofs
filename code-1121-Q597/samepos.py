"""PROVE 11-21: (P2) data. Does the clean journey's rhombus (in M, steps 0..p) ever occupy a position that Jo's rhombus
occupies at some step (as a tile)?  And over the whole clean journey?   usage: python3 samepos.py 5 NW NCH"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def work(w):
    st = Counter(); nstates = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            nstates += 1
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
                    JR = {frozenset(rh(s[1])) for s in Jo} | {frozenset(rh(Jo[-1][2]))}
                    YR = [frozenset(rh(s[1])) for s in A] + [frozenset(rh(A[-1][2]))]
                    st[(tag, 'clean rhombus at a Jo rhombus position, steps 0..p', any(r in JR for r in YR[:p+1]))] += 1
                    st[(tag, 'same, whole clean journey', any(r in JR for r in YR))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns)
    for k in sorted(st, key=str): print(st[k], k)
