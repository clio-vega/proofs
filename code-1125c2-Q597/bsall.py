"""PROVE 2026-11-25 c2: BS(L) on ALL diverging B/C records, by class, plus the MC-target identity.
Per ordered B/C pair (x->B, y->C) at every census state:
  DIV      : which roles diverge (re-run != journey).  Lemma both: one diverges => both do.
  STEP     : both A and A' take a step at p (no early stop), each diverging role.
  L        : roles at a class-1 a-step of their own target at p (Prop orient: exactly one).
  CLASS    : L's re-run type at p:  b -> case b (geometric) ; a_{T'} -> case a ; anything else -> LEFTOVER.
  BS(L)    : bs_test on L's middles.
  DL==DE   : Delta_L == Delta_E, Delta_role = (re-run delta) - (journey delta)  [MC <=> this, given PQ-RIGID]
  COMM     : boards M.Jx.Jy^x == M.Jy.Jx^y
Control: B/B non-commuting pairs: same columns (L via class-1 a-step; BS expected to fail; DL==DE expected to fail).
usage: python3 bsall.py SIZE NW NCH [CHUNKS lo-hi]"""
import os, sys
from collections import Counter
from multiprocessing import Pool
CHARG = sys.argv[4] if len(sys.argv) > 4 else None
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
aT = {'B': 'a1', 'C': 'a0'}
def sub(D, E):
    out = dict(D)
    for k, v in E.items():
        w = tuple(a-b for a, b in zip(out.get(k, (0, 0)), v)); out[k] = w
        if w == (0, 0): del out[k]
    return out
def work(w):
    st = Counter()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            st['states'] += 1
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        comm = tiles_of(b12) == tiles_of(b21)
                        if kind == 'BB' and comm: continue
                        P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                        X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                        Dx, Dy = sub(Xy, X), sub(Yx, Y)
                        roles = [(tx, Jx, Jxy), (ty, Jy, Jyx)]
                        div = [A != Bq for _, A, Bq in roles]
                        tag = kind
                        st[(tag, 'pairs')] += 1
                        if not any(div):
                            st[(tag, 'nondiverging', 'Dx==Dy==0', not Dx and not Dy, 'COMM', comm)] += 1; continue
                        st[(tag, 'DIV both', all(div))] += 1
                        info = []
                        for T, A, Bq in roles:
                            mA, mB, p, q = middles(A, Bq)
                            stepboth = p < len(A) and p < len(Bq)
                            info.append((T, typed(A)[p] if p < len(A) else None, typed(Bq)[p] if p < len(Bq) else None, mA, mB, stepboth))
                        st[(tag, 'STEP both roles', all(i[5] for i in info))] += 1
                        Ls = [r for r in (0, 1) if info[r][1] == (1, aT[info[r][0]])]
                        st[(tag, '#L', len(Ls))] += 1
                        if len(Ls) != 1: continue
                        L = Ls[0]; T, tA, tB, mA, mB, _ = info[L]
                        Tp = roles[1-L][0]
                        if tB == (1, 'b'): cl = 'case b'
                        elif tB == (1, aT[Tp]) and Tp != T: cl = 'case a'
                        else: cl = 'LEFTOVER %s' % (tB,)
                        DL, DE = (Dx, Dy) if L == 0 else (Dy, Dx)
                        st[(tag, cl, 'BS(L)', bs_test(mA, mB), 'DL==DE', DL == DE, 'COMM', comm)] += 1
    return st
if __name__ == '__main__':
    CH = CHARG
    idx = list(range(*map(int, CH.split('-')))) if CH else list(range(NCH))
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for st in pool.imap_unordered(work, idx): tot.update(st)
    print('SIZE', SIZE, 'chunks', CH or 'all', 'of', NCH)
    for k in sorted(tot, key=str): print(tot[k], k)
