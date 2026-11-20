"""PROVE 11-20 c2: out-of-sample check at SIZE 6 (FIXED SUBSAMPLE: chunks 0..K-1 of NCH, census sampler cap), same
journeys as dumpj.py builds (distinct B/C journeys; B/B second journeys).  Checks on real boards:
 (W) Lemma word: every step matches exactly one abstract configuration;
 (AW) ascending R_j-edge translates of tiles removed at j are removed at no other step;
 (F)  ascending translates by ANY unit vector are removed at no later step; control: flat/descending removed later.
usage: python3 aw6.py 6 NW NCH K"""
import sys
K = int(sys.argv[4]); sys.argv = sys.argv[:4]
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
import words2 as W2, words as Wd
from psilib import psi, rh, tr as trf, edgevecs
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
W12 = [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
def jcheck(st, fam, TT, M0, J, D):
    remat = {}
    for k, (_, g, _) in enumerate(J):
        for t in g:
            if t[0] != 'rh': remat.setdefault(t, []).append(k)
    for j, (ty, g, n) in enumerate(J):
        R = rh(g); Ri = W2.it(R); b = min(Ri[1]); ck = Wd.classkey(R)
        remreal = frozenset(W2.it(t) for t in g if t[0] != 'rh'); Rn = W2.it(rh(n))
        hits = [typ for typ in ('a0', 'a1', 'b') if typ[0] == ty and frozenset(W2.tr(t, b) for t in D[TT][(ck, typ)][0]) == remreal
                and min(Rn[1]) == tuple(x+y for x, y in zip(D[TT][(ck, typ)][3], b))]
        st[(fam, '(W) matching configs', len(hits))] += 1
        EV = edgevecs(R)
        for u in g:
            if u[0] == 'rh': continue
            pu = psi(u, TT)
            for v in W12:
                c = 1 if any(v[:2]) else 0
                t = trf(u, v); dp = psi(t, TT)[c] - pu[c]
                lab = {1: 'asc', 0: 'flat', -1: 'desc'}.get(dp, 'other')
                ks = remat.get(t, [])
                later = any(k > j for k in ks); other = any(k != j for k in ks)
                st[(fam, '(F)', lab, 'removed later', later)] += 1
                if lab == 'asc' and v in EV: st[(fam, '(AW) asc R-edge translate removed at another step', other)] += 1
def work(w):
    D = {T: W2.prep(T) for T in 'BC'}
    st = Counter(); seenJ = set(); ns = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1; M0 = bare(b)
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for TT in ('B', 'C'):
                    _, J = mvl(b, x, TT, DV)
                    if not J or tuple(J) in seenJ: continue
                    seenJ.add(tuple(J)); jcheck(st, 'BC', TT, M0, J, D); st['BC journeys'] += 1
                for y in cs:
                    if y == x: continue
                    c1, Kx = mvl(b, x, 'B', DV)
                    if not Kx: continue
                    M1 = apply(M0, Kx)
                    _, Kyx = mvl(c1, y, 'B', DV)
                    if not Kyx or ('BB2', tuple(Kyx)) in seenJ: continue
                    seenJ.add(('BB2', tuple(Kyx))); jcheck(st, 'BB2', 'B', M1, Kyx, D); st['BB2 journeys'] += 1
    st['states'] += ns
    return st
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(K)))
    st = Counter()
    for s in res: st.update(s)
    print('SIZE', SIZE, 'chunks 0..%d of %d' % (K-1, NCH))
    for k in sorted(st, key=str): print(st[k], k)
