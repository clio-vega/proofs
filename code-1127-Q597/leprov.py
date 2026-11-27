"""Provenance of the tiles removed by L^E (L's re-run on M0.J_E) at its middle steps, case b, all k>=0.
If some tile removed by L^E_{m} (m=k+2, the step a+gamma) is created by E_{k+3}, then J_E has the step k+3 (presence).
usage: python3 leprov.py SIZE NW NCH [lo-hi]"""
import sys
src = open("/home/clio/projects/proofs/code-1126-Q597/g3red.py").read()
exec(src.split("\ndef analyse(")[0])
def nonrh(S): return [t for t in S if t[0] != 'rh']
def analyse2(roles, tag, st, M0):
    info = []
    for T, A, Bq in roles:
        mA, mB, p, q = middles(A, Bq)
        info.append(dict(T=T, A=A, Ar=Bq, p=p, la=len(mA), lb=len(mB), tA=typed(A), tB=typed(Bq)))
    Ls = [r for r in (0, 1) if info[r]['p'] < len(info[r]['A']) and info[r]['tA'][info[r]['p']] == (1, aT[info[r]['T']])]
    if len(Ls) != 1: return
    L = Ls[0]; E = 1 - L; iL, iE = info[L], info[E]
    pL, pE, laL = iL['p'], iE['p'], iL['la']
    if iL['tB'][pL] != (1, 'b'): return
    k = laL - 3; kk = 'k=0' if k == 0 else ('k>0' if k > 0 else 'k<0')
    AE, ArL = iE['A'], iL['Ar']
    st[(tag, kk, 'records')] += 1
    st[(tag, kk, 'J_E has step k+3', pE + k + 3 < len(AE))] += 1
    for i in range(0, k + 4):
        if pL + i >= len(ArL): st[(tag, kk, 'L^E short at', i)] += 1; break
        for g in nonrh(ArL[pL + i][1]):
            org = [('E%+d' % (j - pE)) for j, s in enumerate(AE) if g in s[2]]
            if g in M0: org.append('M0')
            st[(tag, kk, 'L^E_%s' % (i if i < k + 1 else 'k+%d' % (i - k)), 'origin', tuple(org))] += 1
def work(w):
    st = Counter()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            M0 = {v[0:1] + (frozenset(v[1]),) for v in b.tiles.values()}
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    tx, ty = 'BC'
                    b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                    b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                    if Jx == Jxy and Jy == Jyx: continue
                    analyse2([(tx, Jx, Jxy), (ty, Jy, Jyx)], 'BC', st, M0)
    return st
if __name__ == '__main__':
    idx = list(range(*map(int, CHARG.split('-')))) if CHARG else list(range(NCH))
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for st in pool.imap_unordered(work, idx): tot.update(st)
    print('SIZE', SIZE, 'chunks', CHARG or 'all', 'of', NCH)
    for k in sorted(tot, key=str): print(tot[k], k)
