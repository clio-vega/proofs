"""PROVE 11-20: out-of-sample check at SIZE 6 (FIXED SUBSAMPLE: chunks 0..K-1 of NCH, cap CAP states per start as in the census sampler).
For every distinct B/C journey and every B/B second journey: REM, CRE, MU at every step (tilepot.check), every candidate in the
abstract table, bad candidates only as undo of k-1.  usage: python3 pot6.py 6 NW NCH K"""
import sys
K0 = int(sys.argv[4]); K = int(sys.argv[5]); sys.argv = sys.argv[:4]
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
import os
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590')
from common import unit_dir
from stepset import candidates
import localenum as LE
from tilepot import check
DVEC = {'B': unit_dir((0, 0, 0, 1)), 'C': unit_dir((0, -1, 0, 0))}
def canon(S, R):
    pts0 = sorted({P[:2] for P in R[1]}); pts1 = sorted({P[2:] for P in R[1]})
    u = tuple(b-a for a, b in zip(*pts0)); base = pts0[0] + pts1[0]
    return u, frozenset((t[0], frozenset(tuple(int(a-b) for a, b in zip(P, base)) for P in t[1])) for t in S)
TAB = {}
for T in 'BC':
    for u in LE.UNITS:
        Rr, Hs = LE.classify(T, u)
        for typ, H in Hs: TAB[(T, u, H)] = LE.analyse(T, Rr, typ, set(H))
def journey_stats(st, fam, TT, M0, J):
    Bd = set(M0); d = DVEC[TT]; prev = None; ok = [True, True, True]
    for k, (typ, g, n) in enumerate(J):
        R = next(t for t in g if t[0] == 'rh')
        for (ty, S, c) in candidates(Bd, R, d):
            u, cs = canon(S, R); ent = TAB.get((TT, u, cs))
            if ent is None: st[(fam, 'candidate NOT in table')] += 1; continue
            bad = ent.get('DL') is False
            st[(fam, 'selected' if S == g else 'other', 'BAD' if bad else 'good', 'undo' if S == prev else ('k=0' if k == 0 else 'not undo'))] += 1
        r = check(TT, g, n)
        for i in range(3): ok[i] &= r[i]
        Bd -= g; Bd |= n; prev = n
    st[(fam, 'journey REM,CRE,MU all steps', tuple(ok))] += 1
def work(w):
    st = Counter(); seen = set(); ns = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1; M0 = bare(b); cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for TT in 'BC':
                    _, J = mvl(b, x, TT, DV)
                    if J and ('BC', tuple(J)) not in seen:
                        seen.add(('BC', tuple(J))); journey_stats(st, 'BC', TT, M0, J)
                for y in cs:
                    if y == x: continue
                    c1, Kx = mvl(b, x, 'B', DV)
                    if not Kx: continue
                    _, Kyx = mvl(c1, y, 'B', DV)
                    if Kyx and ('BB2', tuple(Kyx)) not in seen:
                        seen.add(('BB2', tuple(Kyx))); journey_stats(st, 'BB2', 'B', apply(M0, Kx), Kyx)
    st[('states',)] += ns
    return st
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(K0, K)))
    st = Counter()
    for s in res: st.update(s)
    print("SIZE", SIZE, "chunks %d..%d of %d" % (K0, K-1, NCH))
    for k in sorted(st, key=str): print(st[k], k)
