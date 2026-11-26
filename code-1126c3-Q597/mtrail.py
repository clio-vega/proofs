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
    # alpha as a vertex translation: E^L_i rhombus vs E_i rhombus
    for i in range(2, k + 2):
        s, r = AE[pE+i], ArE[pE+i]
        rs = [t for t in s[1] if t[0] == 'rh']; rr = [t for t in r[1] if t[0] == 'rh']
        if len(rs) != 1 or len(rr) != 1: st[(tag, 'rh count odd')] += 1; continue
        # rhombus created by previous E^L step?
        st[(tag, 'P3 rh of E^L_i created by E^L_{i-1}', rr[0] in ArE[pE+i-1][2])] += 1
        for g in nonrh(r[1]):
            org = []
            for j, sl in enumerate(AL):
                if g in sl[2]: org.append('CL%+d' % (j - pL - i))
            for j in range(pE + i):
                if g in ArE[j][2]: org.append('CEL%+d' % (j - pE - i))
            if g in M0: org.append('M0')
            remL = [j - pL for j, sl in enumerate(AL) if g in sl[1]]
            st[(tag, 'P1 origin(rel i)', tuple(sorted(org)), g[0])] += 1
            st[(tag, 'P2 removed by J_L', bool(remL))] += 1
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
