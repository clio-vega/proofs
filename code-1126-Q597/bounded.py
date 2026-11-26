"""After cancelling the corridor: front/back bounded remainders of Delta_E - Delta_L."""
import pickle, glob
from collections import Counter
from lookstep import delta
st = Counter()
def acc(D, s, sg):
    for k, v in delta(s).items():
        w = tuple(a + sg*b for a, b in zip(D.get(k, (0, 0)), v))
        if w == (0, 0): D.pop(k, None)
        else: D[k] = w
for f in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(f, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        if case == 'a': continue
        k = laL - 3
        F, Bk, FULL = {}, {}, {}
        # front: E^L r0,r1 ; E r0,r1 ; L^E r0 ; L r0,r1
        for i in (0, 1): acc(F, ArE[pE+i], +1); acc(F, AE[pE+i], -1); acc(F, AL[pL+i], +1)
        acc(F, ArL[pL], -1)
        for i in (k+2, k+3): acc(Bk, ArE[pE+i], +1); acc(Bk, AE[pE+i], -1)
        for i in (k+1, k+2): acc(Bk, ArL[pL+i], -1)
        acc(Bk, AL[pL+k+2], +1)
        for i in range(laE): acc(FULL, ArE[pE+i], +1); acc(FULL, AE[pE+i], -1)
        for i in range(laL): acc(FULL, ArL[pL+i], -1); acc(FULL, AL[pL+i], +1)
        tot = dict(F)
        for kk, v in Bk.items():
            w = tuple(a+b for a, b in zip(tot.get(kk, (0, 0)), v))
            if w == (0, 0): tot.pop(kk, None)
            else: tot[kk] = w
        st[('k>0' if k else 'k=0', 'DeltaE-DeltaL==0', not FULL, 'front==0', not F, 'back==0', not Bk, 'front+back==0', not tot, '#front keys', len(F))] += 1
for k in sorted(st, key=str): print(st[k], k)
