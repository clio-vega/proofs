"""Finite part of Lemma win (PROVE 2026-11-26): build the window table (one representative per class up to one Z^4
translation) and check, ON THE REPRESENTATIVES ONLY, the identities (F), (BK) and the case-a identity. Writes windows.pkl.
The identities are linear in delta, and delta(s+v) = T_v delta(s) (key shifted by v's component, displacement unchanged),
so an identity holding on a representative holds on every translate."""
import pickle, glob, sys
from collections import Counter
sys.argv = ['x', '5', '1', '24']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1126-Q597')
from lookstep import delta
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
W = {'front': {}, 'back': {}, 'a': {}}
for fn in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(fn, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        u = [t for t in AE[pE][1] if t[0] == 'sq' and t in AL[pL][1]]; o = min(u[0][1])
        words = lambda J, i, n: ' '.join('%d%s' % t for t in typed(J)[i:i+n])
        if case == 'b':
            k = laL - 3
            fr = tuple(shift(s, o) for s in (AE[pE], AE[pE+1], ArE[pE], ArE[pE+1], AL[pL], AL[pL+1], ArL[pL]))
            W['front'].setdefault((TL, fr), [0, (words(AE, pE, 2), words(ArE, pE, 2), words(AL, pL, 2), words(ArL, pL, 1))])[0] += 1
            ob = min(p for t in AE[pE+k+2][1] for p in t[1])
            bk = tuple(shift(s, ob) for s in (AE[pE+k+2], AE[pE+k+3], ArE[pE+k+2], ArE[pE+k+3], AL[pL+k+2], ArL[pL+k+1], ArL[pL+k+2]))
            W['back'].setdefault((TL, bk), [0, (words(AE, pE+k+2, 2), words(ArE, pE+k+2, 2), words(AL, pL+k+2, 1), words(ArL, pL+k+1, 2))])[0] += 1
        else:
            wa = (tuple(shift(s, o) for s in AE[pE:pE+3]), tuple(shift(s, o) for s in ArE[pE:pE+3]),
                  tuple(shift(s, o) for s in AL[pL:pL+2]), tuple(shift(s, o) for s in ArL[pL:pL+2]))
            W['a'].setdefault((TL, wa), [0, (words(AE, pE, 3), words(ArE, pE, 3), words(AL, pL, 2), words(ArL, pL, 2))])[0] += 1
pickle.dump(W, open('/home/clio/projects/proofs/code-1126-Q597/windows.pkl', 'wb'))
print('FRONT windows (F):  dE^L0+dE^L1-dE0-dE1 - (dL^E0 - dL0 - dL1) = 0')
for (TL, (e0, e1, f0, f1, l0, l1, m0)), (n, wd) in W['front'].items():
    D = comb([(1, f0), (1, f1), (-1, e0), (-1, e1), (-1, m0), (1, l0), (1, l1)])
    print('  T_L=%s  records %3d  E: %s  E^L: %s  L: %s  L^E: %s   F==0: %s' % (TL, n, *wd, not D))
print('BACK windows (BK):  dE^L_{k+2}+dE^L_{k+3}-dE_{k+2}-dE_{k+3} - (dL^E_{k+1}+dL^E_{k+2}-dL_{k+2}) = 0')
for (TL, (e0, e1, f0, f1, l2, m1, m2)), (n, wd) in W['back'].items():
    D = comb([(1, f0), (1, f1), (-1, e0), (-1, e1), (-1, m1), (-1, m2), (1, l2)])
    print('  T_L=%s  records %3d  E: %s  E^L: %s  L: %s  L^E: %s   BK==0: %s' % (TL, n, *wd, not D))
print('CASE a windows:  Delta_E - Delta_L over the whole middles = 0  (and Delta_L = 0, Delta_E = 0 separately)')
for (TL, (Em, ELm, Lm, LEm)), (n, wd) in W['a'].items():
    DE = comb([(1, s) for s in ELm] + [(-1, s) for s in Em]); DL = comb([(1, s) for s in LEm] + [(-1, s) for s in Lm])
    print('  T_L=%s  records %3d  E: %s  E^L: %s  L: %s  L^E: %s   Delta_E==0: %s  Delta_L==0: %s' % (TL, n, *wd, not DE, not DL))
