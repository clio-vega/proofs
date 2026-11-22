"""PROVE 11-22 c2: dump every diverging ordered B/C pair (both roles) and every non-commuting B/B pair (both roles)
with raw step lists, so that the sign analysis runs offline.  Built from code-1122-Q597/pairdich.py, same loop.
Per pair: kind, sideI, I (keys), and per role: (tag, A, B=re-run, Jo) step lists, plus M0 (bare board).
usage: python3 rolesdump.py SIZE NW NCH  (KEEP=k env for the size-6 subsample)"""
import os, pickle
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def work(w):
    out = []; nstates = 0
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
                        P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                        X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                        K = set(X) | set(Y) | set(Yx) | set(Xy)
                        I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                        if kind == 'BB':
                            if tiles_of(b12) == tiles_of(b21): continue
                            roles = (('BBy', Jy, Jyx, Jx), ('BBx', Jx, Jxy, Jy))
                        else:
                            roles = (('x->B', Jx, Jxy, Jy), ('y->C', Jy, Jyx, Jx))
                        divs = [len(middles(A, Bq)[0]) >= 1 and len(middles(A, Bq)[1]) >= 1 for _, A, Bq, _ in roles]
                        if not any(divs): continue
                        out.append(dict(kind=kind, key=key, x=x, y=y, I=I, X=X, Y=Y, Yx=Yx, Xy=Xy, M0=M0,
                                        roles=[(t, A, Bq, Jo) for t, A, Bq, Jo in roles], divs=divs))
    return out, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    allp = []; ns = 0
    for o, n in res: allp += o; ns += n
    pickle.dump(allp, open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'wb'))
    print('SIZE', SIZE, 'states', ns, 'pairs dumped', len(allp))
