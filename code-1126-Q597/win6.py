"""Out-of-sample test of the window classification (W): is every diverging B/C pair's window (front/back in case b,
whole middle in case a) a translate of a class in the size-5 table windows.pkl (FROZEN)?  New classes: check the identity.
usage: python3 win6.py SIZE NW NCH lo-hi"""
import sys, pickle
CHARG = sys.argv[4]
exec(open("/home/clio/projects/proofs/code-1125c2-Q597/bsall.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1126-Q597')
from lookstep import delta
W5 = pickle.load(open('/home/clio/projects/proofs/code-1126-Q597/windows.pkl', 'rb'))
KF = set(W5['front']); KB = set(W5['back']); KA = set(W5['a'])
def shift(s, o):
    f = lambda S: frozenset((t[0], frozenset(tuple(a-b for a, b in zip(p, o)) for p in t[1])) for t in S)
    return (s[0], f(s[1]), f(s[2]))
def comb(terms):
    D = {}
    for sg, s in terms:
        for k, v in delta(s).items():
            w = tuple(a + sg*b for a, b in zip(D.get(k, (0, 0)), v))
            if w == (0, 0): D.pop(k, None)
            else: D[k] = w
    return D
def work(w):
    st = Counter()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    b1, Jx = mvl(b, x, 'B', DV); b2, Jy = mvl(b, y, 'C', DV)
                    b12, Jyx = mvl(b1, y, 'C', DV); b21, Jxy = mvl(b2, x, 'B', DV)
                    if Jx == Jxy and Jy == Jyx: continue
                    info = []
                    for T, A, Bq in (('B', Jx, Jxy), ('C', Jy, Jyx)):
                        mA, mB, p, q = middles(A, Bq); info.append((T, A, Bq, p, len(mA), len(mB)))
                    Ls = [r for r in (0, 1) if typed(info[r][1])[info[r][3]] == (1, aT[info[r][0]])]
                    if len(Ls) != 1: st['#L!=1'] += 1; continue
                    L = Ls[0]; E = 1 - L
                    TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
                    case = 'b' if typed(ArL)[pL] == (1, 'b') else 'a'
                    u = [t for t in AE[pE][1] if t[0] == 'sq' and t in AL[pL][1]]
                    if not u: st[(case, 'no shared square')] += 1; continue
                    o = min(u[0][1])
                    if case == 'b':
                        k = laL - 3
                        if not (laE == laL + 1 and lbE == laE and lbL == laL): st[(case, 'LEN fails')] += 1; continue
                        fr = (TL, tuple(shift(s, o) for s in (AE[pE], AE[pE+1], ArE[pE], ArE[pE+1], AL[pL], AL[pL+1], ArL[pL])))
                        ob = min(p for t in AE[pE+k+2][1] for p in t[1])
                        bk = (TL, tuple(shift(s, ob) for s in (AE[pE+k+2], AE[pE+k+3], ArE[pE+k+2], ArE[pE+k+3], AL[pL+k+2], ArL[pL+k+1], ArL[pL+k+2])))
                        st[(case, 'front in size-5 table', fr in KF)] += 1
                        st[(case, 'back in size-5 table', bk in KB)] += 1
                    else:
                        wa = (TL, (tuple(shift(s, o) for s in AE[pE:pE+laE]), tuple(shift(s, o) for s in ArE[pE:pE+lbE]),
                                   tuple(shift(s, o) for s in AL[pL:pL+laL]), tuple(shift(s, o) for s in ArL[pL:pL+lbL])))
                        st[(case, 'middle window in size-5 table', wa in KA)] += 1
    return st
if __name__ == '__main__':
    idx = list(range(*map(int, CHARG.split('-'))))
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for st in pool.imap_unordered(work, idx): tot.update(st)
    print('SIZE', SIZE, 'chunks', CHARG, 'of', NCH, '(range lo..hi-1)')
    for k in sorted(tot, key=str): print(tot[k], k)
