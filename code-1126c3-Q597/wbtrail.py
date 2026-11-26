"""Merge' trail: origin of each tile removed by E^L_i, 2<=i<=k+1, case b.  usage: python3 mtrail.py SIZE NW NCH [lo-hi]"""
import sys
src = open("/home/clio/projects/proofs/code-1126-Q597/g3red.py").read()
exec(src.split("\ndef analyse(")[0])
def nonrh(S): return [t for t in S if t[0] != 'rh']
def trans(t, v): return (t[0], frozenset(tuple(a + b for a, b in zip(P, v)) for P in t[1]))
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
    k = laL - 3
    AE, ArE, AL = iE['A'], iE['Ar'], iL['A']
    if k <= 0: st[(tag, 'k<=0')] += 1; return
    tL = iL['A']; M = laL
    if pE + k + 3 >= len(ArE) or pE + k + 3 >= len(AE): st[(tag, 'short')] += 1; return
    TE = iE['T']; aTT = aT
    st[(tag, 'Q1 E^L_{k+2}', iE['tB'][pE+k+2], 'T_L', iL['T'])] += 1
    st[(tag, 'Q3 pair', 'E', iE['tA'][pE+k+3], 'EL', iE['tB'][pE+k+3], 'T_L', iL['T'])] += 1
    r = ArE[pE+k+3]; CLk2 = set(AL[pL+k+2][2])
    for g in nonrh(r[1]):
        if g in CLk2: continue
        org = []
        for j, sl in enumerate(AL):
            if g in sl[2]: org.append('CL%+d' % (j - pL))
        if g in M0: org.append('M0')
        remL = any(g in sl[1] for sl in AL)
        st[(tag, 'Q2 E^L_{k+3} tile not in C(L_{k+2}): origin', tuple(org), 'removed by J_L', remL)] += 1
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
