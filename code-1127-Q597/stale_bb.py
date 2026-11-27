"""(Stale_E): at k* = p_E+k+4 (after the back window), D = X* \\ Y* where X* = J_E's board, Y* = E^L's board.
Classify the tiles of D that are NOT created by J_E (those J_E never removes: Lemma fin): origin and the step that took
them off Y.  Then check whether J_E's suffix (steps >= k*) ever removes one (must be no for LEN).  Also LEN itself.
usage: python3 stale_e.py SIZE NW NCH [lo-hi]"""
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
    k = laL - 3; kk = 'k=0' if k == 0 else 'k>0'
    if k < 0: st[(tag, 'k<0')] += 1; return
    AE, ArE, AL = iE['A'], iE['Ar'], iL['A']
    st[(tag, kk, 'LEN', iE['la'] == k + 4 and iE['lb'] == k + 4)] += 1
    ks = pE + k + 4
    if ks > len(AE) or ks > len(ArE): st[(tag, kk, 'short')] += 1; return
    st[(tag, kk, 'W_B pair', 'T_L', iL['T'], 'E', tuple(iE['tA'][pE+k+2:pE+k+4]), 'EL', tuple(iE['tB'][pE+k+2:pE+k+4]))] += 1
    CLk2 = set(AL[pL+k+2][2])
    for g in nonrh(ArE[pE+k+3][1]):
        if g in CLk2: continue
        st[(tag, kk, 'Q2 tile of E^L_{k+3} outside C(L_{k+2}) in M0', g in M0, 'in G(E_{k+3})', g in AE[pE+k+3][1])] += 1
    for g in nonrh(AE[pE+k+3][1]):
        st[(tag, kk, "Q2'' tile of E_{k+3}: removed by E^L at", tuple(j - pE for j in range(len(ArE)) if g in ArE[j][1]),
            'by J_L', any(g in s[1] for s in AL))] += 1
    # removal / creation sets
    remE_before = set().union(*[s[1] for s in AE[:ks]]) if ks else set()
    WL = set().union(*[s[1] for s in AL]) & M0
    remEL = {}
    for j in range(pE, min(ks, len(ArE))):
        for g in ArE[j][1]: remEL.setdefault(g, j - pE)
    creE = set().union(*[s[2] for s in AE])
    remE_after = {}
    for j in range(ks, len(AE)):
        for g in AE[j][1]: remE_after.setdefault(g, j - pE)
    cand = set()
    for g in M0:
        if g in remE_before: continue          # not on X*
        offY = (g in WL) or (g in remEL)
        if not offY: continue
        if g in creE: continue
        cand.add(g)
        kind = ('WL' if g in WL else '') + ('+EL%d' % remEL[g] if g in remEL else '')
        if g in WL:
            js = [j - pL for j, s in enumerate(AL) if g in s[1]]
            kind += '@L%+d' % js[0]
        st[(tag, kk, 'D kind', kind, g[0], 'read by J_E later', g in remE_after)] += 1
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
                    tx, ty = 'BB'
                    if not x < y: continue
                    b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                    b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                    if Jx == Jxy and Jy == Jyx: continue
                    analyse2([(tx, Jx, Jxy), (ty, Jy, Jyx)], 'BB', st, M0)
    return st
if __name__ == '__main__':
    idx = list(range(*map(int, CHARG.split('-')))) if CHARG else list(range(NCH))
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for st in pool.imap_unordered(work, idx): tot.update(st)
    print('SIZE', SIZE, 'chunks', CHARG or 'all', 'of', NCH)
    for k in sorted(tot, key=str): print(tot[k], k)
