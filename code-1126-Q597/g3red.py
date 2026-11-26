"""PROVE 2026-11-26: the corridor reduction of (G3'), checked per diverging pair, B/C and the B/B control.
Per pair (L = role at (1,a_T) at its contact, E the other):
  LEN   : |mid_E| == |mid_L| + 1                                   (no role-swapping symmetry)
  MERGE : J_E[pE+i] == J_L[pL+i] (identical tile sets) exactly for i = 2..k+1, k := |mid_L| - 3 (case b), none in case a
  MERGE': E^L[pE+i] == L^E[pL+i-1] exactly for i = 2..k+1, and == J_E[pE+i] translated by -alpha (alpha from BS(L) vec)
  F, BK : front / back bounded identities (see tex, Lemma red)
  G3'   : Delta_E == Delta_L
  LESION: replace the re-run corridor shift by any other unit vector v != -alpha in the cancellation: does F+BK+corridor still give 0?
usage: python3 g3red.py SIZE NW NCH [lo-hi]"""
import sys
CHARG = sys.argv[4] if len(sys.argv) > 4 else None
exec(open("/home/clio/projects/proofs/code-1125c2-Q597/bsall.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1126-Q597')
from lookstep import delta
from match import canon, vsub
def KEY(s): return (s[0], s[1], s[2])
def acc(D, s, sg, shift=None):
    dd = delta(s)
    for k, v in dd.items():
        if shift is not None:   # translate the key triangle by the side-k component of shift
            off = shift[:2] if k[0] == 0 else shift[2:]
            k = (k[0], frozenset(tuple(a+b for a, b in zip(P, off)) for P in k[1]))
        w = tuple(a + sg*b for a, b in zip(D.get(k, (0, 0)), v))
        if w == (0, 0): D.pop(k, None)
        else: D[k] = w
UNITS = [tuple(int(i == j) * sg for j in range(4)) for i in range(4) for sg in (1, -1)] + \
        [tuple(sg*c for c in u) for u in ((1, -1, 0, 0), (0, 0, 1, -1)) for sg in (1, -1)]
def analyse(roles, tag, st):
    info = []
    for T, A, Bq in roles:
        mA, mB, p, q = middles(A, Bq)
        info.append(dict(T=T, A=A, Ar=Bq, p=p, la=len(mA), lb=len(mB), tA=typed(A), tB=typed(Bq)))
    Ls = [r for r in (0, 1) if info[r]['p'] < len(info[r]['A']) and info[r]['tA'][info[r]['p']] == (1, aT[info[r]['T']])]
    if len(Ls) != 1: st[(tag, '#L', len(Ls))] += 1; return
    L = Ls[0]; E = 1 - L; iL, iE = info[L], info[E]
    pL, pE, laL, laE = iL['p'], iE['p'], iL['la'], iE['la']
    case = 'b' if iL['tB'][pL] == (1, 'b') else 'a'
    k = laL - (3 if case == 'b' else 2)
    st[(tag, case, 'LEN |mid_E|==|mid_L|+1', laE == laL + 1 and iE['lb'] == laE)] += 1
    idxL = {KEY(s): i for i, s in enumerate(iL['A'])}
    sh0 = tuple((i - pE, idxL[KEY(s)] - pL) for i, s in enumerate(iE['A']) if KEY(s) in idxL)
    idxLr = {KEY(s): i for i, s in enumerate(iL['Ar'])}
    sh1 = tuple((i - pE, idxLr[KEY(s)] - pL) for i, s in enumerate(iE['Ar']) if KEY(s) in idxLr)
    want0 = tuple((i, i) for i in range(2, k + 2)); want1 = tuple((i, i - 1) for i in range(2, k + 2))
    st[(tag, case, 'MERGE', sh0 == want0)] += 1
    st[(tag, case, "MERGE'", sh1 == want1)] += 1
    # alpha: L^E[pL] is L[pL+1] translated by -alpha (BS(L)), case b
    malpha = None
    if case == 'b' and laL >= 2:
        a1, m1 = canon(iL['A'][pL+1]); b1, n1 = canon(iL['Ar'][pL])
        if a1 == b1: malpha = vsub(n1, m1)
    if case == 'b' and k > 0 and malpha is not None:
        ok = all(canon(iE['Ar'][pE+i])[0] == canon(iE['A'][pE+i])[0] and vsub(canon(iE['Ar'][pE+i])[1], canon(iE['A'][pE+i])[1]) == malpha for i in range(2, k+2))
        st[(tag, case, "MERGE' corridor shift == -alpha of BS(L)", ok)] += 1
    # Delta_E - Delta_L, front, back
    AE, ArE, AL, ArL = iE['A'], iE['Ar'], iL['A'], iL['Ar']
    FULL = {}
    for i in range(iE['lb']): acc(FULL, ArE[pE+i], +1)
    for i in range(laE): acc(FULL, AE[pE+i], -1)
    for i in range(iL['la']): acc(FULL, AL[pL+i], +1)
    for i in range(iL['lb']): acc(FULL, ArL[pL+i], -1)
    st[(tag, case, "G3' Delta_E==Delta_L", not FULL)] += 1
    if case == 'b' and laE == laL + 1 and iE['lb'] == laE and iL['lb'] == laL:
        F, BK = {}, {}
        for i in (0, 1): acc(F, ArE[pE+i], +1); acc(F, AE[pE+i], -1); acc(F, AL[pL+i], +1)
        acc(F, ArL[pL], -1)
        for i in (k+2, k+3): acc(BK, ArE[pE+i], +1); acc(BK, AE[pE+i], -1)
        for i in (k+1, k+2): acc(BK, ArL[pL+i], -1)
        acc(BK, AL[pL+k+2], +1)
        st[(tag, case, 'F==0', not F, 'BK==0', not BK)] += 1
        # lesion: rebuild E^L's corridor as J_E's corridor shifted by v (instead of the true re-run steps), for each unit v
        if k > 0 and malpha is not None:
            for v in UNITS:
                D = {}
                for i in (0, 1, k+2, k+3): acc(D, ArE[pE+i], +1)
                for i in range(2, k+2): acc(D, AE[pE+i], +1, shift=v)
                for i in range(laE): acc(D, AE[pE+i], -1)  # E's lb==la here
                for i in range(iL['la']): acc(D, AL[pL+i], +1)
                for i in range(iL['lb']): acc(D, ArL[pL+i], -1)
                st[(tag, case, 'LESION corridor shift', 'v==alpha_shift' if v == malpha else 'v!=alpha_shift', 'G3 identity holds', not D)] += 1
                if v != malpha and not D:
                    st[(tag, case, 'LESION pass with v!=alpha: v - alpha has zero component on the corridor key side', 'alpha', malpha, 'v', v)] += 1
def work(w):
    st = Counter()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        if kind == 'BB' and tiles_of(b12) == tiles_of(b21): continue
                        if Jx == Jxy and Jy == Jyx: continue
                        analyse([(tx, Jx, Jxy), (ty, Jy, Jyx)], kind, st)
    return st
if __name__ == '__main__':
    idx = list(range(*map(int, CHARG.split('-')))) if CHARG else list(range(NCH))
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for st in pool.imap_unordered(work, idx): tot.update(st)
    print('SIZE', SIZE, 'chunks', CHARG or 'all', 'of', NCH)
    for k in sorted(tot, key=str): print(tot[k], k)
