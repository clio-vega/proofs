"""Window anatomy of (SW) + the 'triangle' test (T): tau;x->B;b->C == tau;x->C ?   (size-5 closure, first NSEED seeds)
usage: python3 swwin.py NSEED"""
import sys
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
def mvn(b, X, ij, dvec):
    LOG.clear(); IBC.STEPLOG.clear(); t = b.nest[X]['pos'].pop(ij); R, Y, ij2, s = b.journey(t, X, unit_dir(dvec)); return (Y, ij2), list(LOG)
st = Counter(); NSEED = int(sys.argv[1]); ex = []
for k, (key, b0) in enumerate(census()):
    if k >= NSEED: break
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            # triangle test on every corner
            bB = copy.deepcopy(b); (YB, bcell), K = mvn(bB, 0, x, dirs['B'])
            bC = copy.deepcopy(b); (YC, ccell), H = mvn(bC, 0, x, dirs['C'])
            Bn = bB.nest[1]
            for nm, dv in (('wB', Bn['w']), ('uB', Bn['u'])):
                t = copy.deepcopy(bB)
                try: (Y2, c2), J = mvn(t, 1, bcell, dv)
                except Exception as e: st['T '+nm+' err'] += 1; continue
                st[f'T {nm} lands nest {Y2}'] += 1
                if Y2 == 2: st[f'T {nm} ->C equals tau;x->C'] += tiles_of(t) == tiles_of(bC)
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            L = copy.deepcopy(bB); _, G = mvn(L, 0, y, dirs['C'])
            R = copy.deepcopy(bC); _, Gp = mvn(R, 0, y, dirs['B'])
            st['SW n'] += 1; st['SW same'] += tiles_of(L) == tiles_of(R)
            sC = suffix_len(G, H); sB = suffix_len(Gp, K)
            iC = len(G) - sC; iB = len(Gp) - sB
            st[('C-side pre-merge, skipped', iC, len(H) - sC)] += 1
            st[('B-side pre-merge, skipped', iB, len(K) - sB)] += 1
            # left window = K[:len(K)-sB] + G[:iC] ; right window = H[:len(H)-sC] + Gp[:iB]
            WL = K[:len(K)-sB] + G[:iC]; WR = H[:len(H)-sC] + Gp[:iB]
            vL = set(v for s_ in WL for vs in s_ for v in vs); vR = set(v for s_ in WR for vs in s_ for v in vs)
            st['window vertex sets equal'] += vL == vR
            st[('window lengths L,R', len(WL), len(WR))] += 1
for k in sorted(st, key=str): print(st[k], k)
