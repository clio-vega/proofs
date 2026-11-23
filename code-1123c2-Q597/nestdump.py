"""PROVE 11-23 c2: rolesdump.py loop (unchanged), B/C diverging pairs only, plus the nest data of M0 (X, u, w, pos per nest),
so that the engine's termination test (Board.packed_at + is_ideal) can be evaluated offline at any rhombus.
usage: python3 nestdump.py SIZE NW NCH  (KEEP=k env for the size-6 subsample)"""
import os, pickle
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def work(w):
    out = []; nstates = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            nstates += 1; M0 = bare(b)
            cs = corners(set(b.nest[0]['pos']))
            nests = [dict(X=info['X'], u=info['u'], w=info['w'], pos=dict(info['pos'])) for info in b.nest]
            for x in cs:
                for y in cs:
                    if x == y: continue
                    b1, Jx = mvl(b, x, 'B', DV); b2, Jy = mvl(b, y, 'C', DV)
                    b12, Jyx = mvl(b1, y, 'C', DV); b21, Jxy = mvl(b2, x, 'B', DV)
                    roles = (('x->B', Jx, Jxy, Jy), ('y->C', Jy, Jyx, Jx))
                    divs = [len(middles(A, Bq)[0]) >= 1 and len(middles(A, Bq)[1]) >= 1 for _, A, Bq, _ in roles]
                    if not any(divs): continue
                    nb = lambda bb: [dict(pos=dict(info['pos'])) for info in bb.nest]
                    out.append(dict(key=key, x=x, y=y, M0=M0, nests=nests, nests1=nb(b1), nests2=nb(b2),
                                    commute=tiles_of(b12) == tiles_of(b21),
                                    roles=[(t, A, Bq, Jo) for t, A, Bq, Jo in roles], divs=divs))
    return out, nstates
if __name__ == '__main__':
    idx = list(range(int(os.environ.get('KEEP', NCH))))
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, idx))
    allp = []; ns = 0
    for o, n in res: allp += o; ns += n
    pickle.dump(allp, open('/home/clio/projects/proofs/code-1123c2-Q597/nest%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'wb'))
    print('SIZE', SIZE, 'states', ns, 'pairs dumped', len(allp))
