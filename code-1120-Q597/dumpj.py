"""PROVE 11-20: dump every distinct B/C journey (<=5 census, same sampler as 11-19 c2), plus B/B second journeys
(journey of y->B after x->B, x!=y corners) as the control.  usage: python3 dumpj.py 5 NW NCH  -> dumpj5-<w>.pkl"""
import sys, pickle
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
def work(w):
    out = {'BC': [], 'BB2': []}; seenJ = set(); ns = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1; M0 = bare(b)
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for TT in ('B', 'C'):
                    _, J = mvl(b, x, TT, DV)
                    if not J or tuple(J) in seenJ: continue
                    seenJ.add(tuple(J)); out['BC'].append((TT, frozenset(M0), list(J)))
                for y in cs:
                    if y == x: continue
                    c1, Kx = mvl(b, x, 'B', DV)
                    if not Kx: continue
                    M1 = apply(M0, Kx)
                    _, Kyx = mvl(c1, y, 'B', DV)
                    if not Kyx or ('BB2', tuple(Kyx)) in seenJ: continue
                    seenJ.add(('BB2', tuple(Kyx))); out['BB2'].append(('B', frozenset(M1), list(Kyx)))
    pickle.dump(out, open('dumpj%d-%d.pkl' % (SIZE, w), 'wb'))
    return ns, len(out['BC']), len(out['BB2'])
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    print('SIZE', SIZE, 'states', sum(r[0] for r in res), 'BC journeys', sum(r[1] for r in res), 'BB2', sum(r[2] for r in res))
